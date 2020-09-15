//
//  vdb_apic_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/13.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "vdb_apic_solver3.hpp"
#include "vdb_point_searcher3.hpp"
#include "vdb_samplers3.h"

using namespace vdb;

ApicSolver3::ApicSolver3()
: ApicSolver3({1, 1, 1}, {1, 1, 1}, {0, 0, 0}) {
}

ApicSolver3::ApicSolver3(const vox::Size3& resolution,
                         const vox::Vector3D& gridSpacing,
                         const vox::Vector3D& gridOrigin)
: PicSolver3(resolution, gridSpacing, gridOrigin) {
}

ApicSolver3::~ApicSolver3(){}

void ApicSolver3::transferFromParticlesToGrids() {
    auto flow = vdbSystemData()->velocity();
    const auto particles = particleSystemData();
    const auto positions = particles->positions();
    auto velocities = particles->velocities();
    const size_t numberOfParticles = particles->numberOfParticles();
    const auto hh = flow->gridSpacing() / 2.0;
    const auto bbox = flow->boundingBox();
    
    // Allocate buffers
    _cX.resize(numberOfParticles);
    _cY.resize(numberOfParticles);
    _cZ.resize(numberOfParticles);
    
    // Clear velocity to zero
    flow->clear();
    
    const auto uPos = vdbSystemData()->velocity()->uPosition();
    const auto vPos = vdbSystemData()->velocity()->vPosition();
    const auto wPos = vdbSystemData()->velocity()->wPosition();
    openvdb::DoubleGrid::Ptr uWeight = openvdb::DoubleGrid::create(0.0);
    openvdb::DoubleGrid::Ptr vWeight = openvdb::DoubleGrid::create(0.0);
    openvdb::DoubleGrid::Ptr wWeight = openvdb::DoubleGrid::create(0.0);
    _uMarkers.resize(flow->uSize());
    _vMarkers.resize(flow->vSize());
    _wMarkers.resize(flow->wSize());
    _uMarkers.set(0);
    _vMarkers.set(0);
    _wMarkers.set(0);
    LinearGridSampler<openvdb::DoubleGrid> uSampler(vdbSystemData()->velocity()->getUGrid(),
                                                    flow->uSize(),
                                                    flow->gridSpacing(),
                                                    flow->uOrigin());
    LinearGridSampler<openvdb::DoubleGrid> vSampler(vdbSystemData()->velocity()->getVGrid(),
                                                    flow->vSize(),
                                                    flow->gridSpacing(),
                                                    flow->vOrigin());
    LinearGridSampler<openvdb::DoubleGrid> wSampler(vdbSystemData()->velocity()->getWGrid(),
                                                    flow->wSize(),
                                                    flow->gridSpacing(),
                                                    flow->wOrigin());
    
    for (size_t i = 0; i < numberOfParticles; ++i) {
        std::array<openvdb::Coord, 8> indices;
        std::array<double, 8> weights;
        
        auto uPosClamped = positions[i];
        uPosClamped.y = vox::clamp(
                                   uPosClamped.y,
                                   bbox.lowerCorner.y + hh.y,
                                   bbox.upperCorner.y - hh.y);
        uPosClamped.z = vox::clamp(
                                   uPosClamped.z,
                                   bbox.lowerCorner.z + hh.z,
                                   bbox.upperCorner.z - hh.z);
        
        uSampler.getCoordinatesAndWeights(uPosClamped, &indices, &weights);
        for (int j = 0; j < 8; ++j) {
            vox::Vector3D gridPos = uPos(indices[j]);
            double apicTerm = _cX[i].dot(gridPos - uPosClamped);
            
            vdbSystemData()->velocity()
            ->getUGrid()->tree().setValueOn(indices[j], vdbSystemData()->velocity()
                                            ->getUGrid()->tree().getValue(indices[j])
                                            + weights[j] * (velocities[i].x + apicTerm));
            uWeight->tree().setValueOn(indices[j], uWeight->tree().getValue(indices[j]) + weights[j]);
            _uMarkers(indices[j].x(), indices[j].y(), indices[j].z()) = 1;
        }
        
        auto vPosClamped = positions[i];
        vPosClamped.x = vox::clamp(
                                   vPosClamped.x,
                                   bbox.lowerCorner.x + hh.x,
                                   bbox.upperCorner.x - hh.x);
        vPosClamped.z = vox::clamp(
                                   vPosClamped.z,
                                   bbox.lowerCorner.z + hh.z,
                                   bbox.upperCorner.z - hh.z);
        
        vSampler.getCoordinatesAndWeights(vPosClamped, &indices, &weights);
        for (int j = 0; j < 8; ++j) {
            vox::Vector3D gridPos = vPos(indices[j]);
            double apicTerm = _cY[i].dot(gridPos - vPosClamped);
            
            vdbSystemData()->velocity()
            ->getVGrid()->tree().setValueOn(indices[j], vdbSystemData()->velocity()
                                            ->getVGrid()->tree().getValue(indices[j])
                                            + weights[j] * (velocities[i].y + apicTerm));
            vWeight->tree().setValueOn(indices[j], vWeight->tree().getValue(indices[j]) + weights[j]);
            _vMarkers(indices[j].x(), indices[j].y(), indices[j].z()) = 1;
        }
        
        auto wPosClamped = positions[i];
        wPosClamped.x = vox::clamp(
                                   wPosClamped.x,
                                   bbox.lowerCorner.x + hh.x,
                                   bbox.upperCorner.x - hh.x);
        wPosClamped.y = vox::clamp(
                                   wPosClamped.y,
                                   bbox.lowerCorner.y + hh.y,
                                   bbox.upperCorner.y - hh.y);
        
        wSampler.getCoordinatesAndWeights(wPosClamped, &indices, &weights);
        for (int j = 0; j < 8; ++j) {
            vox::Vector3D gridPos = wPos(indices[j]);
            double apicTerm = _cZ[i].dot(gridPos - wPosClamped);
            
            vdbSystemData()->velocity()
            ->getWGrid()->tree().setValueOn(indices[j], vdbSystemData()->velocity()
                                            ->getWGrid()->tree().getValue(indices[j])
                                            + weights[j] * (velocities[i].z + apicTerm));
            wWeight->tree().setValueOn(indices[j], wWeight->tree().getValue(indices[j]) + weights[j]);
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

void ApicSolver3::transferFromGridsToParticles() {
    const auto flow = vdbSystemData()->velocity();
    const auto particles = particleSystemData();
    auto positions = particles->positions();
    auto velocities = particles->velocities();
    const size_t numberOfParticles = particles->numberOfParticles();
    const auto hh = vdbSystemData()->velocity()->gridSpacing() / 2.0;
    const auto bbox = flow->boundingBox();
    
    // Allocate buffers
    _cX.resize(numberOfParticles);
    _cY.resize(numberOfParticles);
    _cZ.resize(numberOfParticles);
    _cX.set(vox::Vector3D());
    _cY.set(vox::Vector3D());
    _cZ.set(vox::Vector3D());
    LinearGridSampler<openvdb::DoubleGrid> uSampler(vdbSystemData()->velocity()->getUGrid(),
                                                    flow->uSize(),
                                                    flow->gridSpacing(),
                                                    flow->uOrigin());
    LinearGridSampler<openvdb::DoubleGrid> vSampler(vdbSystemData()->velocity()->getVGrid(),
                                                    flow->vSize(),
                                                    flow->gridSpacing(),
                                                    flow->vOrigin());
    LinearGridSampler<openvdb::DoubleGrid> wSampler(vdbSystemData()->velocity()->getWGrid(),
                                                    flow->wSize(),
                                                    flow->gridSpacing(),
                                                    flow->wOrigin());
    
    vox::parallelFor(vox::kZeroSize, numberOfParticles, [&](size_t i) {
        velocities[i] = vdbSystemData()->velocity()->sample(positions[i]);
        
        std::array<openvdb::Coord, 8> indices;
        std::array<vox::Vector3D, 8> gradWeights;
        
        // x
        auto uPosClamped = positions[i];
        uPosClamped.y = vox::clamp(
                                   uPosClamped.y,
                                   bbox.lowerCorner.y + hh.y,
                                   bbox.upperCorner.y - hh.y);
        uPosClamped.z = vox::clamp(
                                   uPosClamped.z,
                                   bbox.lowerCorner.z + hh.z,
                                   bbox.upperCorner.z - hh.z);
        uSampler.getCoordinatesAndGradientWeights(uPosClamped,
                                                  &indices, &gradWeights);
        for (int j = 0; j < 8; ++j) {
            _cX[i] += gradWeights[j] * vdbSystemData()->velocity()
            ->getUGrid()->tree().getValue(indices[j]);
        }
        
        // y
        auto vPosClamped = positions[i];
        vPosClamped.x = vox::clamp(
                                   vPosClamped.x,
                                   bbox.lowerCorner.x + hh.x,
                                   bbox.upperCorner.x - hh.x);
        vPosClamped.z = vox::clamp(
                                   vPosClamped.z,
                                   bbox.lowerCorner.z + hh.z,
                                   bbox.upperCorner.z - hh.z);
        vSampler.getCoordinatesAndGradientWeights(vPosClamped,
                                                  &indices, &gradWeights);
        for (int j = 0; j < 8; ++j) {
            _cY[i] += gradWeights[j] * vdbSystemData()->velocity()
            ->getVGrid()->tree().getValue(indices[j]);
        }
        
        // z
        auto wPosClamped = positions[i];
        wPosClamped.x = vox::clamp(
                                   wPosClamped.x,
                                   bbox.lowerCorner.x + hh.x,
                                   bbox.upperCorner.x - hh.x);
        wPosClamped.y = vox::clamp(
                                   wPosClamped.y,
                                   bbox.lowerCorner.y + hh.y,
                                   bbox.upperCorner.y - hh.y);
        wSampler.getCoordinatesAndGradientWeights(wPosClamped,
                                                  &indices, &gradWeights);
        for (int j = 0; j < 8; ++j) {
            _cZ[i] += gradWeights[j] * vdbSystemData()->velocity()
            ->getWGrid()->tree().getValue(indices[j]);
        }
    });
}

//-----------------------------------------------
ApicSolver3::Builder ApicSolver3::builder() {
    return Builder();
}

ApicSolver3 ApicSolver3::Builder::build() const {
    return ApicSolver3(_resolution, getGridSpacing(), _gridOrigin);
}

ApicSolver3Ptr ApicSolver3::Builder::makeShared() const {
    return std::shared_ptr<ApicSolver3>(
                                        new ApicSolver3(_resolution,
                                                        getGridSpacing(),
                                                        _gridOrigin),
                                        [] (ApicSolver3* obj) {
        delete obj;
    });
}
