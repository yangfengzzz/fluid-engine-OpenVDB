//
//  vdb_flip_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/13.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "vdb_flip_solver3.hpp"
#include "vdb_samplers3.h"
#include "../src.common/pch.h"

using namespace vdb;

FlipSolver3::FlipSolver3() : FlipSolver3({1, 1, 1}, {1, 1, 1}, {0, 0, 0}) {}

FlipSolver3::FlipSolver3(const vox::Size3& resolution,
                         const vox::Vector3D& gridSpacing,
                         const vox::Vector3D& gridOrigin)
: PicSolver3(resolution, gridSpacing, gridOrigin) {
}

FlipSolver3::~FlipSolver3() {}

double FlipSolver3::picBlendingFactor() const { return _picBlendingFactor; }

void FlipSolver3::setPicBlendingFactor(double factor) {
    _picBlendingFactor = vox::clamp(factor, 0.0, 1.0);
}

void FlipSolver3::transferFromParticlesToGrids() {
    PicSolver3::transferFromParticlesToGrids();
    
    _uDelta = openvdb::DoubleGrid::create(0.0);
    _uDelta->topologyUnion(*vdbSystemData()->velocity()->getUGrid());
    _uDelta->setTransform(vdbSystemData()->velocity()->getUGrid()->transform().copy());
    
    _vDelta = openvdb::DoubleGrid::create(0.0);
    _vDelta->topologyUnion(*vdbSystemData()->velocity()->getVGrid());
    _vDelta->setTransform(vdbSystemData()->velocity()->getVGrid()->transform().copy());
    
    _wDelta = openvdb::DoubleGrid::create(0.0);
    _wDelta->topologyUnion(*vdbSystemData()->velocity()->getWGrid());
    _wDelta->setTransform(vdbSystemData()->velocity()->getWGrid()->transform().copy());
    
    vdbSystemData()->velocity()->parallelForEachUIndex([&](const openvdb::Coord& coord) {
        _uDelta->tree().setValue(coord, vdbSystemData()->velocity()
                                 ->getUGrid()->tree().getValue(coord));
    });
    vdbSystemData()->velocity()->parallelForEachVIndex([&](const openvdb::Coord& coord) {
        _vDelta->tree().setValue(coord, vdbSystemData()->velocity()
                                 ->getVGrid()->tree().getValue(coord));
    });
    vdbSystemData()->velocity()->parallelForEachWIndex([&](const openvdb::Coord& coord) {
        _wDelta->tree().setValue(coord, vdbSystemData()->velocity()
                                 ->getWGrid()->tree().getValue(coord));
    });
}

void FlipSolver3::transferFromGridsToParticles() {
    auto flow = vdbSystemData()->velocity();
    auto positions = particleSystemData()->positions();
    auto velocities = particleSystemData()->velocities();
    size_t numberOfParticles = particleSystemData()->numberOfParticles();
    
    // Compute delta
    _uDelta->topologyUnion(*vdbSystemData()->velocity()->getUGrid());
    vdbSystemData()->velocity()->parallelForEachUIndex([&](const openvdb::Coord& coord) {
        _uDelta->tree().setValueOn(coord, vdbSystemData()->velocity()
                                   ->getUGrid()->tree().getValue(coord)
                                   - _uDelta->tree().getValue(coord) );
    });
    
    _vDelta->topologyUnion(*vdbSystemData()->velocity()->getVGrid());
    vdbSystemData()->velocity()->parallelForEachVIndex([&](const openvdb::Coord& coord) {
        _vDelta->tree().setValueOn(coord, vdbSystemData()->velocity()
                                   ->getVGrid()->tree().getValue(coord)
                                   - _vDelta->tree().getValue(coord) );
    });
    
    _wDelta->topologyUnion(*vdbSystemData()->velocity()->getWGrid());
    vdbSystemData()->velocity()->parallelForEachWIndex([&](const openvdb::Coord& coord) {
        _wDelta->tree().setValueOn(coord, vdbSystemData()->velocity()
                                   ->getWGrid()->tree().getValue(coord)
                                   - _wDelta->tree().getValue(coord) );
    });
    
    LinearGridSampler<openvdb::DoubleGrid> uSampler(_uDelta,
                                                    flow->uSize(),
                                                    flow->gridSpacing(),
                                                    flow->uOrigin());
    LinearGridSampler<openvdb::DoubleGrid> vSampler(_vDelta, flow->vSize(),
                                                    flow->gridSpacing(),
                                                    flow->vOrigin());
    LinearGridSampler<openvdb::DoubleGrid> wSampler(_wDelta, flow->wSize(),
                                                    flow->gridSpacing(),
                                                    flow->wOrigin());
    auto sampler = [uSampler, vSampler, wSampler](const vox::Vector3D& x) {
        double u = uSampler(x);
        double v = vSampler(x);
        double w = wSampler(x);
        return vox::Vector3D(u, v, w);
    };
    
    // Transfer delta to the particles
    vox::parallelFor(vox::kZeroSize, numberOfParticles, [&](size_t i) {
        vox::Vector3D flipVel = velocities[i] + sampler(positions[i]);
        if (_picBlendingFactor > 0.0) {
            vox::Vector3D picVel = vdbSystemData()->velocity()->sample(positions[i]);
            flipVel = lerp(flipVel, picVel, _picBlendingFactor);
        }
        velocities[i] = flipVel;
    });
}

//-----------------------------------------------
FlipSolver3::Builder FlipSolver3::builder() {
    return Builder();
}

FlipSolver3 FlipSolver3::Builder::build() const {
    return FlipSolver3(_resolution, getGridSpacing(), _gridOrigin);
}

FlipSolver3Ptr FlipSolver3::Builder::makeShared() const {
    return std::shared_ptr<FlipSolver3>(
                                        new FlipSolver3(_resolution,
                                                        getGridSpacing(),
                                                        _gridOrigin),
                                        [] (FlipSolver3* obj) {
        delete obj;
    });
}
