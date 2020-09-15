//
//  vdb_pic_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/6.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_samplers3.h"
#include "vdb_helper.h"
#include "../src.common/level_set_utils.h"
#include "vdb_pic_solver3.hpp"
#include "vdb_point_searcher3.hpp"
#include "../src.common/timer.h"

#include <algorithm>

using namespace vdb;

PicSolver3::PicSolver3()
: PicSolver3({1, 1, 1}, {1, 1, 1}, {0, 0, 0}) {
}

PicSolver3::PicSolver3(
                       const vox::Size3& resolution,
                       const vox::Vector3D& gridSpacing,
                       const vox::Vector3D& gridOrigin)
: FluidSolver3(resolution, gridSpacing, gridOrigin) {
    auto grids = vdbSystemData();
    _signedDistanceFieldId = grids->addScalarData(
                                                  std::make_shared<CellCenteredScalarGrid3::Builder>(), vox::kMaxD);
    _particles = std::make_shared<ParticleSystemData3>();
}

PicSolver3::~PicSolver3() {
}

ScalarGrid3Ptr PicSolver3::signedDistanceField() const {
    return vdbSystemData()->scalarDataAt(_signedDistanceFieldId);
}

const ParticleSystemData3Ptr& PicSolver3::particleSystemData() const {
    return _particles;
}

const ParticleEmitter3Ptr& PicSolver3::particleEmitter() const {
    return _particleEmitter;
}

void PicSolver3::setParticleEmitter(const ParticleEmitter3Ptr& newEmitter) {
    _particleEmitter = newEmitter;
    newEmitter->setTarget(_particles);
}
//------------------------------------------------------------------------------
void PicSolver3::onInitialize() {
    FluidSolver3::onInitialize();
    
    vox::Timer timer;
    updateParticleEmitter(0.0);
    LOG(INFO) << "Update particle emitter took "
    << timer.durationInSeconds() << " seconds";
}

void PicSolver3::onBeginAdvanceTimeStep(double timeIntervalInSeconds) {
    UNUSED_VARIABLE(timeIntervalInSeconds);
    
    LOG(INFO) << "Number of PIC-type particles: "
    << _particles->numberOfParticles();
    
    vox::Timer timer;
    updateParticleEmitter(timeIntervalInSeconds);
    LOG(INFO) << "Update particle emitter took "
    << timer.durationInSeconds() << " seconds";
    
    LOG(INFO) << "Number of PIC-type particles: "
    << _particles->numberOfParticles();
    
    timer.reset();
    transferFromParticlesToGrids();
    LOG(INFO) << "transferFromParticlesToGrids took "
    << timer.durationInSeconds() << " seconds";
    
    timer.reset();
    buildSignedDistanceField();
    LOG(INFO) << "buildSignedDistanceField took "
    << timer.durationInSeconds() << " seconds";
    
    timer.reset();
    extrapolateVelocityToAir();
    LOG(INFO) << "extrapolateVelocityToAir took "
    << timer.durationInSeconds() << " seconds";
    
    applyBoundaryCondition();
}

void PicSolver3::computeAdvection(double timeIntervalInSeconds) {
    vox::Timer timer;
    extrapolateVelocityToAir();
    LOG(INFO) << "extrapolateVelocityToAir took "
    << timer.durationInSeconds() << " seconds";
    
    applyBoundaryCondition();
    
    timer.reset();
    transferFromGridsToParticles();
    LOG(INFO) << "transferFromGridsToParticles took "
    << timer.durationInSeconds() << " seconds";
    
    timer.reset();
    moveParticles(timeIntervalInSeconds);
    LOG(INFO) << "moveParticles took "
    << timer.durationInSeconds() << " seconds";
}

vox::ScalarField3Ptr PicSolver3::fluidSdf() const {
    return signedDistanceField();
}

//------------------------------------------------------------------------------
void PicSolver3::transferFromParticlesToGrids() {
    auto flow = vdbSystemData()->velocity();
    auto positions = _particles->positions();
    auto velocities = _particles->velocities();
    size_t numberOfParticles = _particles->numberOfParticles();
    
    // Clear velocity to zero
    flow->clear();
    
    // Weighted-average velocity
    openvdb::DoubleGrid::Ptr uWeight = openvdb::DoubleGrid::create(0.0);
    openvdb::DoubleGrid::Ptr vWeight = openvdb::DoubleGrid::create(0.0);
    openvdb::DoubleGrid::Ptr wWeight = openvdb::DoubleGrid::create(0.0);
    LinearGridSampler<openvdb::DoubleGrid> uSampler(flow->getUGrid(),
                                                    flow->uSize(),
                                                    flow->gridSpacing(),
                                                    flow->uOrigin());
    LinearGridSampler<openvdb::DoubleGrid> vSampler(flow->getVGrid(),
                                                    flow->vSize(),
                                                    flow->gridSpacing(),
                                                    flow->vOrigin());
    LinearGridSampler<openvdb::DoubleGrid> wSampler(flow->getWGrid(),
                                                    flow->wSize(),
                                                    flow->gridSpacing(),
                                                    flow->wOrigin());
    _uMarkers.resize(flow->uSize());
    _vMarkers.resize(flow->vSize());
    _wMarkers.resize(flow->wSize());
    _uMarkers.set(0);
    _vMarkers.set(0);
    _wMarkers.set(0);
    
    for (size_t i = 0; i < numberOfParticles; ++i) {
        std::array<openvdb::Coord, 8> indices;
        std::array<double, 8> weights;
        
        uSampler.getCoordinatesAndWeights(positions[i],
                                          &indices, &weights);
        for (int j = 0; j < 8; ++j) {
            vdbSystemData()->velocity()
            ->getUGrid()->tree().setValueOn(indices[j],
                                            vdbSystemData()->velocity()
                                            ->getUGrid()->tree().getValue(indices[j])
                                            + velocities[i].x * weights[j]);
            uWeight->tree().setValueOn(indices[j],
                                       uWeight->tree().getValue(indices[j])
                                       + weights[j]);
            _uMarkers(indices[j].x(), indices[j].y(), indices[j].z()) = 1;
        }
        
        vSampler.getCoordinatesAndWeights(positions[i],
                                          &indices, &weights);
        for (int j = 0; j < 8; ++j) {
            vdbSystemData()->velocity()
            ->getVGrid()->tree().setValueOn(indices[j], vdbSystemData()->velocity()
                                            ->getVGrid()->tree().getValue(indices[j])
                                            + velocities[i].y * weights[j]);
            vWeight->tree().setValueOn(indices[j],
                                       vWeight->tree().getValue(indices[j])
                                       + weights[j]);
            _vMarkers(indices[j].x(), indices[j].y(), indices[j].z()) = 1;
        }
        
        wSampler.getCoordinatesAndWeights(positions[i],
                                          &indices, &weights);
        for (int j = 0; j < 8; ++j) {
            vdbSystemData()->velocity()
            ->getWGrid()->tree().setValueOn(indices[j], vdbSystemData()->velocity()
                                            ->getWGrid()->tree().getValue(indices[j])
                                            + velocities[i].z * weights[j]);
            wWeight->tree().setValueOn(indices[j],
                                       wWeight->tree().getValue(indices[j])
                                       + weights[j]);
            _wMarkers(indices[j].x(), indices[j].y(), indices[j].z()) = 1;
        }
    }
    
    GridHelper<openvdb::DoubleGrid> uHelper(uWeight);
    uHelper.parallelForEachDataPointIndex([&](const openvdb::Coord &coord)->void{
        double val = uWeight->tree().getValue(coord);
        if (val > 0.0) {
            vdbSystemData()->velocity()
            ->getUGrid()->tree().setValueOn(coord, vdbSystemData()->velocity()
                                            ->getUGrid()->tree().getValue(coord)/val);
        }
    });
    
    GridHelper<openvdb::DoubleGrid> vHelper(vWeight);
    vHelper.parallelForEachDataPointIndex([&](const openvdb::Coord &coord)->void{
        double val= vWeight->tree().getValue(coord);
        if (val > 0.0) {
            vdbSystemData()->velocity()
            ->getVGrid()->tree().setValueOn(coord, vdbSystemData()->velocity()
                                            ->getVGrid()->tree().getValue(coord)/val);
        }
    });
    
    GridHelper<openvdb::DoubleGrid> wHelper(wWeight);
    wHelper.parallelForEachDataPointIndex([&](const openvdb::Coord &coord)->void{
        double val = wWeight->tree().getValue(coord);
        if (val > 0.0) {
            vdbSystemData()->velocity()
            ->getWGrid()->tree().setValueOn(coord, vdbSystemData()->velocity()
                                            ->getWGrid()->tree().getValue(coord)/val);
        }
    });
}

void PicSolver3::transferFromGridsToParticles() {
    auto positions = _particles->positions();
    auto velocities = _particles->velocities();
    size_t numberOfParticles = _particles->numberOfParticles();
    
    vox::parallelFor(vox::kZeroSize, numberOfParticles, [&](size_t i) {
        velocities[i] = vdbSystemData()->velocity()->sample(positions[i]);
    });
}

void PicSolver3::moveParticles(double timeIntervalInSeconds) {
    auto flow = vdbSystemData()->velocity();
    auto positions = _particles->positions();
    auto velocities = _particles->velocities();
    size_t numberOfParticles = _particles->numberOfParticles();
    int domainBoundaryFlag = closedDomainBoundaryFlag();
    vox::BoundingBox3D boundingBox = flow->boundingBox();
    
    vox::parallelFor(vox::kZeroSize, numberOfParticles, [&](size_t i) {
        vox::Vector3D pt0 = positions[i];
        vox::Vector3D pt1 = pt0;
        vox::Vector3D vel = velocities[i];
        
        // Adaptive time-stepping
        unsigned int numSubSteps
        = static_cast<unsigned int>(std::max(maxCfl(), 1.0));
        double dt = timeIntervalInSeconds / numSubSteps;
        for (unsigned int t = 0; t < numSubSteps; ++t) {
            vox::Vector3D vel0 = vdbSystemData()->velocity()->sample(pt0);
            
            // Mid-point rule
            vox::Vector3D midPt = pt0 + 0.5 * dt * vel0;
            vox::Vector3D midVel = vdbSystemData()->velocity()->sample(midPt);
            pt1 = pt0 + dt * midVel;
            
            pt0 = pt1;
        }
        
        if ((domainBoundaryFlag & vox::kDirectionLeft)
            && pt1.x <= boundingBox.lowerCorner.x) {
            pt1.x = boundingBox.lowerCorner.x;
            vel.x = 0.0;
        }
        if ((domainBoundaryFlag & vox::kDirectionRight)
            && pt1.x >= boundingBox.upperCorner.x) {
            pt1.x = boundingBox.upperCorner.x;
            vel.x = 0.0;
        }
        if ((domainBoundaryFlag & vox::kDirectionDown)
            && pt1.y <= boundingBox.lowerCorner.y) {
            pt1.y = boundingBox.lowerCorner.y;
            vel.y = 0.0;
        }
        if ((domainBoundaryFlag & vox::kDirectionUp)
            && pt1.y >= boundingBox.upperCorner.y) {
            pt1.y = boundingBox.upperCorner.y;
            vel.y = 0.0;
        }
        if ((domainBoundaryFlag & vox::kDirectionBack)
            && pt1.z <= boundingBox.lowerCorner.z) {
            pt1.z = boundingBox.lowerCorner.z;
            vel.z = 0.0;
        }
        if ((domainBoundaryFlag & vox::kDirectionFront)
            && pt1.z >= boundingBox.upperCorner.z) {
            pt1.z = boundingBox.upperCorner.z;
            vel.z = 0.0;
        }
        
        positions[i] = pt1;
        velocities[i] = vel;
    });
    
    vox::Collider3Ptr col = collider();
    if (col != nullptr) {
        vox::parallelFor(
                         vox::kZeroSize,
                         numberOfParticles,
                         [&](size_t i) {
            openvdb::Vec3d pt(positions[i].x,
                              positions[i].y,
                              positions[i].z);
            openvdb::Vec3d vel(velocities[i].x,
                               velocities[i].y,
                               velocities[i].z);
            col->resolveCollision(
                                  0.0, 0.0,
                                  &positions[i],
                                  &velocities[i]);
        });
    }
}

void PicSolver3::extrapolateVelocityToAir() {
    unsigned int depth = static_cast<unsigned int>(std::ceil(maxCfl()));
    extrapolateToRegion<openvdb::DoubleGrid>(vdbSystemData()->velocity()->getUGrid(),
                                             _uMarkers, depth,
                                             vdbSystemData()->velocity()->getUGrid());
    extrapolateToRegion<openvdb::DoubleGrid>(vdbSystemData()->velocity()->getVGrid(),
                                             _vMarkers, depth,
                                             vdbSystemData()->velocity()->getVGrid());
    extrapolateToRegion<openvdb::DoubleGrid>(vdbSystemData()->velocity()->getWGrid(),
                                             _wMarkers, depth,
                                             vdbSystemData()->velocity()->getWGrid());
}

void PicSolver3::buildSignedDistanceField() {
    auto sdf = signedDistanceField();
    double maxH = vox::max3(
                            sdf->gridSpacing().x,
                            sdf->gridSpacing().y,
                            sdf->gridSpacing().z);
    double radius = 1.2 * maxH / std::sqrt(2.0);
    double sdfBandRadius = 2.0 * radius;
    
    _particles->buildNeighborSearcher(2 * radius);
    auto searcher = _particles->neighborSearcher();
    signedDistanceField()->fill([&] (const vox::Vector3D& pt)->double {
        double minDist = sdfBandRadius;
        searcher->forEachNearbyPoint(pt,
                                     sdfBandRadius, [&] (size_t, const vox::Vector3D& x) {
            minDist = std::min(minDist, pt.distanceTo(x));
        });
        return minDist - radius;
    }, vox::ExecutionPolicy::kSerial);
    
    extrapolateIntoCollider(signedDistanceField().get());
}

void PicSolver3::updateParticleEmitter(double timeIntervalInSeconds) {
    if (_particleEmitter != nullptr) {
        _particleEmitter->update(currentTimeInSeconds(), timeIntervalInSeconds);
    }
}

//-------------------------------------------------------------------------------
PicSolver3::Builder PicSolver3::builder() {
    return Builder();
}

PicSolver3 PicSolver3::Builder::build() const {
    return PicSolver3(_resolution, getGridSpacing(), _gridOrigin);
}

PicSolver3Ptr PicSolver3::Builder::makeShared() const {
    return std::shared_ptr<PicSolver3>(
                                       new PicSolver3(_resolution,
                                                      getGridSpacing(),
                                                      _gridOrigin),
                                       [] (PicSolver3* obj) {
        delete obj;
    });
}
