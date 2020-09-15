//
//  vdb_fluid_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/6.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "../src.common/array_utils.h"
#include "../src.common/constant_scalar_field3.h"
#include "../src.common/constants.h"
#include "vdb_backward_euler_diffusion_solver3.hpp"
#include "vdb_blocked_boundary_condition_solver3.hpp"
#include "vdb_fluid_solver3.hpp"
#include "vdb_single_phase_pressure_solver3.hpp"
#include "vdb_fractional_single_phase_pressure_solver3.hpp"
#include "../src.common/level_set_utils.h"
#include "../src.common/surface_to_implicit3.h"
#include "../src.common/timer.h"

#include <algorithm>

using namespace vdb;

FluidSolver3::FluidSolver3()
: FluidSolver3({1, 1, 1}, {1, 1, 1}, {0, 0, 0}) {}

FluidSolver3::FluidSolver3(const vox::Size3& resolution,
                           const vox::Vector3D& gridSpacing,
                           const vox::Vector3D& gridOrigin) {
    _grids = std::make_shared<GridSystemData3>();
    _grids->resize(resolution, gridSpacing, gridOrigin);
    
    setDiffusionSolver(std::make_shared<BackwardEulerDiffusionSolver3>());
    setPressureSolver(
                      std::make_shared<FractionalSinglePhasePressureSolver3>());
    setIsUsingFixedSubTimeSteps(false);
}

FluidSolver3::~FluidSolver3() {}

const vox::Vector3D& FluidSolver3::gravity() const { return _gravity; }

void FluidSolver3::setGravity(const vox::Vector3D& newGravity) {
    _gravity = newGravity;
}

double FluidSolver3::viscosityCoefficient() const {
    return _viscosityCoefficient;
}

void FluidSolver3::setViscosityCoefficient(double newValue) {
    _viscosityCoefficient = std::max(newValue, 0.0);
}

double FluidSolver3::cfl(double timeIntervalInSeconds) const {
    auto vel = _grids->velocity();
    double maxVel = 0.0;
    vel->forEachCellIndex([&](uint i, uint j, uint k) {
        vox::Vector3D v = vel->valueAtCellCenter(i, j, k)
        + timeIntervalInSeconds * _gravity;
        maxVel = std::max(maxVel, v.x);
        maxVel = std::max(maxVel, v.y);
        maxVel = std::max(maxVel, v.z);
    });
    
    vox::Vector3D gridSpacing = _grids->gridSpacing();
    double minGridSize = vox::min3(gridSpacing.x,
                                   gridSpacing.y,
                                   gridSpacing.z);
    
    return maxVel * timeIntervalInSeconds / minGridSize;
}

double FluidSolver3::maxCfl() const { return _maxCfl; }

void FluidSolver3::setMaxCfl(double newCfl) {
    _maxCfl = std::max(newCfl, vox::kEpsilonD);
}

bool FluidSolver3::useCompressedLinearSystem() const {
    return _useCompressedLinearSys;
}

void FluidSolver3::setUseCompressedLinearSystem(bool onoff) {
    _useCompressedLinearSys = onoff;
}

//-------------------------------------------------------------------------
const AdvectionSolver3Ptr& FluidSolver3::advectionSolver() const {
    return _advectionSolver;
}

void FluidSolver3::setAdvectionSolver(const AdvectionSolver3Ptr& newSolver) {
    _advectionSolver = newSolver;
}

const DiffusionSolver3Ptr& FluidSolver3::diffusionSolver() const {
    return _diffusionSolver;
}

void FluidSolver3::setDiffusionSolver(
                                      const DiffusionSolver3Ptr& newSolver) {
    _diffusionSolver = newSolver;
}

const PressureSolver3Ptr& FluidSolver3::pressureSolver() const {
    return _pressureSolver;
}

void FluidSolver3::setPressureSolver(
                                     const PressureSolver3Ptr& newSolver) {
    _pressureSolver = newSolver;
    if (_pressureSolver != nullptr) {
        _boundaryConditionSolver =
        _pressureSolver->suggestedBoundaryConditionSolver();
        
        // Apply domain boundary flag
        _boundaryConditionSolver->setClosedDomainBoundaryFlag(
                                                              _closedDomainBoundaryFlag);
    }
}

int FluidSolver3::closedDomainBoundaryFlag() const {
    return _closedDomainBoundaryFlag;
}

void FluidSolver3::setClosedDomainBoundaryFlag(int flag) {
    _closedDomainBoundaryFlag = flag;
    _boundaryConditionSolver->setClosedDomainBoundaryFlag(
                                                          _closedDomainBoundaryFlag);
}

const GridSystemData3Ptr& FluidSolver3::vdbSystemData() const {
    return _grids;
}

void FluidSolver3::resizeGrid(const vox::Size3& newSize,
                              const vox::Vector3D& newGridSpacing,
                              const vox::Vector3D& newGridOrigin) {
    _grids->resize(newSize, newGridSpacing, newGridOrigin);
}

vox::Size3 FluidSolver3::resolution() const { return _grids->resolution(); }

vox::Vector3D FluidSolver3::gridSpacing() const {
    return _grids->gridSpacing(); }

vox::Vector3D FluidSolver3::gridOrigin() const { return _grids->origin(); }

const FaceCenteredGrid3Ptr& FluidSolver3::velocity() const {
    return _grids->velocity();
}

const vox::Collider3Ptr& FluidSolver3::collider() const { return _collider; }

void FluidSolver3::setCollider(const vox::Collider3Ptr& newCollider) {
    _collider = newCollider;
}

const vox::GridEmitter3Ptr& FluidSolver3::emitter() const { return _emitter; }

void FluidSolver3::setEmitter(const vox::GridEmitter3Ptr& newEmitter) {
    _emitter = newEmitter;
}

//----------------------------------------------------------------
void FluidSolver3::onInitialize() {
    // When initializing the solver, update the collider and emitter state as
    // well since they also affects the initial condition of the simulation.
    vox::Timer timer;
    updateCollider(0.0);
    LOG(INFO) << "Update collider took " << timer.durationInSeconds()
    << " seconds";
    
    timer.reset();
    updateEmitter(0.0);
    LOG(INFO) << "Update emitter took " << timer.durationInSeconds()
    << " seconds";
}

void FluidSolver3::onAdvanceTimeStep(double timeIntervalInSeconds) {
    // The minimum grid resolution is 1x1.
    if (_grids->resolution().x == 0 || _grids->resolution().y == 0 ||
        _grids->resolution().z == 0) {
        LOG(WARNING) << "Empty grid. Skipping the simulation.";
        return;
    }
    
    beginAdvanceTimeStep(timeIntervalInSeconds);
    
    vox::Timer timer;
    computeExternalForces(timeIntervalInSeconds);
    LOG(INFO) << "Computing external force took " << timer.durationInSeconds()
    << " seconds";
    
    timer.reset();
    computeViscosity(timeIntervalInSeconds);
    LOG(INFO) << "Computing viscosity force took " << timer.durationInSeconds()
    << " seconds";
    
    timer.reset();
    computePressure(timeIntervalInSeconds);
    LOG(INFO) << "Computing pressure force took " << timer.durationInSeconds()
    << " seconds";
    
    timer.reset();
    computeAdvection(timeIntervalInSeconds);
    LOG(INFO) << "Computing advection force took " << timer.durationInSeconds()
    << " seconds";
    
    endAdvanceTimeStep(timeIntervalInSeconds);
}

unsigned int FluidSolver3::numberOfSubTimeSteps(
                                                double timeIntervalInSeconds) const {
    double currentCfl = cfl(timeIntervalInSeconds);
    return static_cast<unsigned int>(
                                     std::max(std::ceil(currentCfl / _maxCfl), 1.0));
}

void FluidSolver3::onBeginAdvanceTimeStep(double timeIntervalInSeconds) {
    UNUSED_VARIABLE(timeIntervalInSeconds);
}

void FluidSolver3::onEndAdvanceTimeStep(double timeIntervalInSeconds) {
    UNUSED_VARIABLE(timeIntervalInSeconds);
}

void FluidSolver3::computeExternalForces(double timeIntervalInSeconds) {
    computeGravity(timeIntervalInSeconds);
}

void FluidSolver3::computeViscosity(double timeIntervalInSeconds) {
    if (_diffusionSolver != nullptr && _viscosityCoefficient > vox::kEpsilonD) {
        auto vel = velocity();
        auto vel0 = std::dynamic_pointer_cast<FaceCenteredGrid3>(vel->clone());
        
        _diffusionSolver->solve(*vel0, _viscosityCoefficient,
                                timeIntervalInSeconds, vel.get(),
                                *colliderSdf(), *fluidSdf());
        applyBoundaryCondition();
    }
}

void FluidSolver3::computePressure(double timeIntervalInSeconds) {
    if (_pressureSolver != nullptr) {
        auto vel = velocity();
        auto vel0 = std::dynamic_pointer_cast<FaceCenteredGrid3>(vel->clone());
        
        _pressureSolver->solve(*vel0, timeIntervalInSeconds, vel.get(),
                               *colliderSdf(), *colliderVelocityField(),
                               *fluidSdf(), _useCompressedLinearSys);
        applyBoundaryCondition();
    }
}

void FluidSolver3::computeAdvection(double timeIntervalInSeconds) {
    UNUSED_VARIABLE(timeIntervalInSeconds);
}

void FluidSolver3::advect(double timeIntervalInSeconds) {
    auto vel = velocity();
    if (_advectionSolver != nullptr) {
        // Solve advections for custom scalar fields
        size_t n = _grids->numberOfAdvectableScalarData();
        for (size_t i = 0; i < n; ++i) {
            auto grid = _grids->advectableScalarDataAt(i);
            auto grid0 = grid->clone();
            _advectionSolver->advect(*grid0, *vel, timeIntervalInSeconds,
                                     grid.get(), *colliderSdf());
            extrapolateIntoCollider(grid.get());
        }
        
        // Solve advections for custom vector fields
        n = _grids->numberOfAdvectableVectorData();
        size_t velIdx = _grids->velocityIndex();
        for (size_t i = 0; i < n; ++i) {
            // Handle velocity layer separately.
            if (i == velIdx) {
                continue;
            }
            
            auto grid = _grids->advectableVectorDataAt(i);
            auto grid0 = grid->clone();
            
            auto collocated =
            std::dynamic_pointer_cast<CollocatedVectorGrid3>(grid);
            auto collocated0 =
            std::dynamic_pointer_cast<CollocatedVectorGrid3>(grid0);
            if (collocated != nullptr) {
                _advectionSolver->advect(*collocated0, *vel,
                                         timeIntervalInSeconds,
                                         collocated.get(), *colliderSdf());
                extrapolateIntoCollider(collocated.get());
                continue;
            }
        }
        applyBoundaryCondition();
    }
}

void FluidSolver3::computeGravity(double timeIntervalInSeconds) {
    if (_gravity.lengthSquared() > vox::kEpsilonD) {
        if (std::abs(_gravity.x) > vox::kEpsilonD) {
            _grids->velocity()->forEachUIndex([&](const openvdb::Coord& coord) {
                _grids->velocity()
                ->getUGrid()->tree().setValue(coord, _grids->velocity()
                                              ->getUGrid()->tree().getValue(coord)
                                              + timeIntervalInSeconds * _gravity.x );
            });
        }
        
        if (std::abs(_gravity.y) > vox::kEpsilonD) {
            _grids->velocity()->forEachVIndex([&](const openvdb::Coord& coord) {
                _grids->velocity()
                ->getVGrid()->tree().setValue(coord, _grids->velocity()
                                              ->getVGrid()->tree().getValue(coord)
                                              + timeIntervalInSeconds * _gravity.y );
            });
        }
        
        if (std::abs(_gravity.z) > vox::kEpsilonD) {
            _grids->velocity()->forEachWIndex([&](const openvdb::Coord& coord) {
                _grids->velocity()
                ->getWGrid()->tree().setValue(coord, _grids->velocity()
                                              ->getWGrid()->tree().getValue(coord)
                                              + timeIntervalInSeconds * _gravity.z );
            });
        }
        
        applyBoundaryCondition();
    }
}

void FluidSolver3::applyBoundaryCondition() {
    auto vel = _grids->velocity();
    
    if (vel != nullptr && _boundaryConditionSolver != nullptr) {
        unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
        _boundaryConditionSolver->constrainVelocity(vel.get(), depth);
    }
}

void FluidSolver3::extrapolateIntoCollider(ScalarGrid3* grid) {
    vox::Array3<char> marker(grid->dataSize());
    auto pos = grid->dataPosition();
    marker.parallelForEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(colliderSdf()->sample(pos(openvdb::Coord(i, j, k))))) {
            marker(i, j, k) = 0;
        } else {
            marker(i, j, k) = 1;
        }
    });
    
    unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
    extrapolateToRegion<openvdb::DoubleGrid>(grid->getGrid(),
                                             marker, depth,
                                             grid->getGrid());
}

void FluidSolver3::extrapolateIntoCollider(CollocatedVectorGrid3* grid) {
    vox::Array3<char> marker(grid->dataSize());
    auto pos = grid->dataPosition();
    marker.parallelForEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(colliderSdf()->sample(pos(openvdb::Coord(i, j, k))))) {
            marker(i, j, k) = 0;
        } else {
            marker(i, j, k) = 1;
        }
    });
    
    unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
    extrapolateToRegion<openvdb::Vec3dGrid>(grid->getGrid(),
                                            marker, depth,
                                            grid->getGrid());
}

void FluidSolver3::extrapolateIntoCollider(FaceCenteredGrid3* grid) {
    auto uPos = grid->uPosition();
    auto vPos = grid->vPosition();
    auto wPos = grid->wPosition();
    
    vox::Array3<char> uMarker(grid->uSize() );
    vox::Array3<char> vMarker(grid->vSize() );
    vox::Array3<char> wMarker(grid->wSize() );
    
    uMarker.parallelForEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(colliderSdf()->sample(uPos(openvdb::Coord(i, j, k))))) {
            uMarker(i, j, k) = 0;
        } else {
            uMarker(i, j, k) = 1;
        }
    });
    
    vMarker.parallelForEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(colliderSdf()->sample(vPos(openvdb::Coord(i, j, k))))) {
            vMarker(i, j, k) = 0;
        } else {
            vMarker(i, j, k) = 1;
        }
    });
    
    wMarker.parallelForEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(colliderSdf()->sample(wPos(openvdb::Coord(i, j, k))))) {
            wMarker(i, j, k) = 0;
        } else {
            wMarker(i, j, k) = 1;
        }
    });
    
    unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
    extrapolateToRegion<openvdb::DoubleGrid>(grid->getUGrid(),
                                             uMarker, depth,
                                             grid->getUGrid());
    extrapolateToRegion<openvdb::DoubleGrid>(grid->getVGrid(),
                                             vMarker, depth,
                                             grid->getVGrid());
    extrapolateToRegion<openvdb::DoubleGrid>(grid->getWGrid(),
                                             wMarker, depth,
                                             grid->getWGrid());
}

vox::ScalarField3Ptr FluidSolver3::fluidSdf() const {
    return std::make_shared<vox::ConstantScalarField3>(-vox::kMaxD);
}

vox::ScalarField3Ptr FluidSolver3::colliderSdf() const {
    return _boundaryConditionSolver->colliderSdf();
}

vox::VectorField3Ptr FluidSolver3::colliderVelocityField() const {
    return _boundaryConditionSolver->colliderVelocityField();
}

void FluidSolver3::beginAdvanceTimeStep(double timeIntervalInSeconds) {
    // Update collider and emitter
    vox::Timer timer;
    updateCollider(timeIntervalInSeconds);
    LOG(INFO) << "Update collider took " << timer.durationInSeconds()
    << " seconds";
    
    timer.reset();
    updateEmitter(timeIntervalInSeconds);
    LOG(INFO) << "Update emitter took " << timer.durationInSeconds()
    << " seconds";
    
    // Update boundary condition solver
    if (_boundaryConditionSolver != nullptr) {
        _boundaryConditionSolver->updateCollider(
                                                 _collider,
                                                 _grids->resolution(),
                                                 _grids->gridSpacing(),
                                                 _grids->origin());
    }
    
    // Apply boundary condition to the velocity field in case the field got
    // updated externally.
    applyBoundaryCondition();
    
    // Invoke callback
    onBeginAdvanceTimeStep(timeIntervalInSeconds);
}

void FluidSolver3::endAdvanceTimeStep(double timeIntervalInSeconds) {
    // Invoke callback
    onEndAdvanceTimeStep(timeIntervalInSeconds);
}

void FluidSolver3::updateCollider(double timeIntervalInSeconds) {
    if (_collider != nullptr) {
        _collider->update(currentTimeInSeconds(), timeIntervalInSeconds);
    }
}

void FluidSolver3::updateEmitter(double timeIntervalInSeconds) {
    if (_emitter != nullptr) {
        _emitter->update(currentTimeInSeconds(), timeIntervalInSeconds);
    }
}

FluidSolver3::Builder FluidSolver3::builder() { return Builder(); }
