//
//  vdb_particle_system_data3.cpp
//  src.vdb
//
//  Created by Feng Yang on 2020/5/19.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "../src.common/pch.h"

#include "../src.common/parallel.h"
#include "vdb_particle_system_data3.hpp"
#include "vdb_point_searcher3.hpp"
#include "../src.common/timer.h"

#include <algorithm>
#include <vector>

using namespace vdb;

ParticleSystemData3::ParticleSystemData3()
: ParticleSystemData3(0) {
}

ParticleSystemData3::ParticleSystemData3(size_t numberOfParticles) {
    _positionIdx = addVectorData();
    _velocityIdx = addVectorData();
    _forceIdx = addVectorData();
    
    // Use PointParallelHashGridSearcher3 by default
    _neighborSearcher = std::make_shared<PointVDBSearcher3>();
    
    resize(numberOfParticles);
}

ParticleSystemData3::ParticleSystemData3(const ParticleSystemData3& other) {
    set(other);
}

ParticleSystemData3::~ParticleSystemData3() {
}

void ParticleSystemData3::resize(size_t newNumberOfParticles) {
    _numberOfParticles = newNumberOfParticles;
    
    for (auto& attr : _scalarDataList) {
        attr.resize(newNumberOfParticles, 0.0);
    }
    
    for (auto& attr : _vectorDataList) {
        attr.resize(newNumberOfParticles, vox::Vector3D());
    }
}

size_t ParticleSystemData3::numberOfParticles() const {
    return _numberOfParticles;
}

size_t ParticleSystemData3::addScalarData(double initialVal) {
    size_t attrIdx = _scalarDataList.size();
    _scalarDataList.emplace_back(numberOfParticles(), initialVal);
    return attrIdx;
}

size_t ParticleSystemData3::addVectorData(const vox::Vector3D& initialVal) {
    size_t attrIdx = _vectorDataList.size();
    _vectorDataList.emplace_back(numberOfParticles(), initialVal);
    return attrIdx;
}

double ParticleSystemData3::radius() const {
    return _radius;
}

void ParticleSystemData3::setRadius(double newRadius) {
    _radius = std::max(newRadius, 0.0);
}

double ParticleSystemData3::mass() const {
    return _mass;
}

void ParticleSystemData3::setMass(double newMass) {
    _mass = std::max(newMass, 0.0);
}

vox::ConstArrayAccessor1<vox::Vector3D> ParticleSystemData3::positions() const {
    return vectorDataAt(_positionIdx);
}

vox::ArrayAccessor1<vox::Vector3D> ParticleSystemData3::positions() {
    return vectorDataAt(_positionIdx);
}

vox::ConstArrayAccessor1<vox::Vector3D> ParticleSystemData3::velocities() const {
    return vectorDataAt(_velocityIdx);
}

vox::ArrayAccessor1<vox::Vector3D> ParticleSystemData3::velocities() {
    return vectorDataAt(_velocityIdx);
}

vox::ConstArrayAccessor1<vox::Vector3D> ParticleSystemData3::forces() const {
    return vectorDataAt(_forceIdx);
}

vox::ArrayAccessor1<vox::Vector3D> ParticleSystemData3::forces() {
    return vectorDataAt(_forceIdx);
}

vox::ConstArrayAccessor1<double> ParticleSystemData3::scalarDataAt(
                                                                   size_t idx) const {
    return _scalarDataList[idx].constAccessor();
}

vox::ArrayAccessor1<double> ParticleSystemData3::scalarDataAt(size_t idx) {
    return _scalarDataList[idx].accessor();
}

vox::ConstArrayAccessor1<vox::Vector3D> ParticleSystemData3::vectorDataAt(
                                                                          size_t idx) const {
    return _vectorDataList[idx].constAccessor();
}

vox::ArrayAccessor1<vox::Vector3D> ParticleSystemData3::vectorDataAt(size_t idx) {
    return _vectorDataList[idx].accessor();
}

void ParticleSystemData3::addParticle(
                                      const vox::Vector3D& newPosition,
                                      const vox::Vector3D& newVelocity,
                                      const vox::Vector3D& newForce) {
    vox::Array1<vox::Vector3D> newPositions = {newPosition};
    vox::Array1<vox::Vector3D> newVelocities = {newVelocity};
    vox::Array1<vox::Vector3D> newForces = {newForce};
    
    addParticles(
                 newPositions.constAccessor(),
                 newVelocities.constAccessor(),
                 newForces.constAccessor());
}

void ParticleSystemData3::addParticles(
                                       const vox::ConstArrayAccessor1<vox::Vector3D>& newPositions,
                                       const vox::ConstArrayAccessor1<vox::Vector3D>& newVelocities,
                                       const vox::ConstArrayAccessor1<vox::Vector3D>& newForces) {
    JET_THROW_INVALID_ARG_IF(
                             newVelocities.size() > 0
                             && newVelocities.size() != newPositions.size());
    JET_THROW_INVALID_ARG_IF(
                             newForces.size() > 0 && newForces.size() != newPositions.size());
    
    size_t oldNumberOfParticles = numberOfParticles();
    size_t newNumberOfParticles = oldNumberOfParticles + newPositions.size();
    
    resize(newNumberOfParticles);
    
    auto pos = positions();
    auto vel = velocities();
    auto frc = forces();
    
    vox::parallelFor(vox::kZeroSize, newPositions.size(),
                     [&](size_t i) {
        pos[i + oldNumberOfParticles] = newPositions[i];
    });
    
    if (newVelocities.size() > 0) {
        vox::parallelFor(vox::kZeroSize, newPositions.size(),
                         [&](size_t i) {
            vel[i + oldNumberOfParticles] = newVelocities[i];
        });
    }
    
    if (newForces.size() > 0) {
        vox::parallelFor(vox::kZeroSize, newPositions.size(),
                         [&](size_t i) {
            frc[i + oldNumberOfParticles] = newForces[i];
        });
    }
}

const vox::PointNeighborSearcher3Ptr& ParticleSystemData3::neighborSearcher() const {
    return _neighborSearcher;
}

void ParticleSystemData3::setNeighborSearcher(
                                              const vox::PointNeighborSearcher3Ptr& newNeighborSearcher) {
    _neighborSearcher = newNeighborSearcher;
}

const std::vector<std::vector<size_t>>&
ParticleSystemData3::neighborLists() const {
    return _neighborLists;
}

void ParticleSystemData3::buildNeighborSearcher(double maxSearchRadius) {
    vox::Timer timer;
    
    // Use PointParallelHashGridSearcher3 by default
    _neighborSearcher = std::make_shared<PointVDBSearcher3>();
    
    _neighborSearcher->build(positions());
    
    LOG(INFO) << "Building neighbor searcher took: "
    << timer.durationInSeconds()
    << " seconds";
}

void ParticleSystemData3::buildNeighborLists(double maxSearchRadius) {
    vox::Timer timer;
    
    _neighborLists.resize(numberOfParticles());
    
    auto points = positions();
    for (size_t i = 0; i < numberOfParticles(); ++i) {
        vox::Vector3D origin = points[i];
        _neighborLists[i].clear();
        
        _neighborSearcher->forEachNearbyPoint(
                                              origin,
                                              maxSearchRadius,
                                              [&](size_t j, const vox::Vector3D&) {
            if (i != j) {
                _neighborLists[i].push_back(j);
            }
        });
    }
    
    LOG(INFO) << "Building neighbor list took: "
    << timer.durationInSeconds()
    << " seconds";
}

void ParticleSystemData3::set(const ParticleSystemData3& other) {
    _radius = other._radius;
    _mass = other._mass;
    _positionIdx = other._positionIdx;
    _velocityIdx = other._velocityIdx;
    _forceIdx = other._forceIdx;
    _numberOfParticles = other._numberOfParticles;
    
    for (auto& attr : other._scalarDataList) {
        _scalarDataList.emplace_back(attr);
    }
    
    for (auto& attr : other._vectorDataList) {
        _vectorDataList.emplace_back(attr);
    }
    
    _neighborSearcher = other._neighborSearcher->clone();
    _neighborLists = other._neighborLists;
}

ParticleSystemData3& ParticleSystemData3::operator=(
                                                    const ParticleSystemData3& other) {
    set(other);
    return *this;
}
