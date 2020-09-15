//
//  vdb_smoke_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/21.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "../src.common/serial.h"
#include "vdb_smoke_solver3.hpp"
#include "vdb_semi_lagrangian3.hpp"

#include <algorithm>

using namespace vdb;

GridSmokeSolver3::GridSmokeSolver3()
: FlipSolver3({1, 1, 1}, {1, 1, 1}, {0, 0, 0}) {}

GridSmokeSolver3::GridSmokeSolver3(const vox::Size3& resolution,
                                   const vox::Vector3D& gridSpacing,
                                   const vox::Vector3D& gridOrigin)
: FlipSolver3(resolution, gridSpacing, gridOrigin) {
    auto grids = vdbSystemData();
    _smokeDensityDataId
    = grids->addAdvectableScalarData(
                                     std::make_shared<CellCenteredScalarGrid3::Builder>(), 0.0);
    vox::serialFor(
                   vox::kZeroSize, resolution.x,
                   vox::kZeroSize, resolution.y,
                   vox::kZeroSize, resolution.z,
                   [&](uint i, uint j, uint k) {
        smokeDensity()->getGrid()->tree().setValueOn(openvdb::Coord(i, j, k));
    });
    
    _temperatureDataId
    = grids->addAdvectableScalarData(
                                     std::make_shared<CellCenteredScalarGrid3::Builder>(), 0.0);
    vox::serialFor(
                   vox::kZeroSize, resolution.x,
                   vox::kZeroSize, resolution.y,
                   vox::kZeroSize, resolution.z,
                   [&](uint i, uint j, uint k) {
        temperature()->getGrid()->tree().setValueOn(openvdb::Coord(i, j, k));
    });
    
    setAdvectionSolver(std::make_shared<SemiLagrangian3>());
}

GridSmokeSolver3::~GridSmokeSolver3() {}

double GridSmokeSolver3::smokeDiffusionCoefficient() const {
    return _smokeDiffusionCoefficient;
}

void GridSmokeSolver3::setSmokeDiffusionCoefficient(double newValue) {
    _smokeDiffusionCoefficient = std::max(newValue, 0.0);
}

double GridSmokeSolver3::temperatureDiffusionCoefficient() const {
    return _temperatureDiffusionCoefficient;
}

void GridSmokeSolver3::setTemperatureDiffusionCoefficient(double newValue) {
    _temperatureDiffusionCoefficient = std::max(newValue, 0.0);
}

double GridSmokeSolver3::buoyancySmokeDensityFactor() const {
    return _buoyancySmokeDensityFactor;
}

void GridSmokeSolver3::setBuoyancySmokeDensityFactor(double newValue) {
    _buoyancySmokeDensityFactor = newValue;
}

double GridSmokeSolver3::buoyancyTemperatureFactor() const {
    return _buoyancyTemperatureFactor;
}

void GridSmokeSolver3::setBuoyancyTemperatureFactor(double newValue) {
    _buoyancyTemperatureFactor = newValue;
}

double GridSmokeSolver3::smokeDecayFactor() const { return _smokeDecayFactor; }

void GridSmokeSolver3::setSmokeDecayFactor(double newValue) {
    _smokeDecayFactor = vox::clamp(newValue, 0.0, 1.0);
}

double GridSmokeSolver3::smokeTemperatureDecayFactor() const {
    return _temperatureDecayFactor;
}

void GridSmokeSolver3::setTemperatureDecayFactor(double newValue) {
    _temperatureDecayFactor = vox::clamp(newValue, 0.0, 1.0);
}

ScalarGrid3Ptr GridSmokeSolver3::smokeDensity() const {
    return vdbSystemData()->advectableScalarDataAt(_smokeDensityDataId);
}

ScalarGrid3Ptr GridSmokeSolver3::temperature() const {
    return vdbSystemData()->advectableScalarDataAt(_temperatureDataId);
}

void GridSmokeSolver3::onEndAdvanceTimeStep(double timeIntervalInSeconds) {
    computeDiffusion(timeIntervalInSeconds);
}

void GridSmokeSolver3::computeExternalForces(double timeIntervalInSeconds) {
    computeBuoyancyForce(timeIntervalInSeconds);
}

void GridSmokeSolver3::computeAdvection(double timeIntervalInSeconds){
    FlipSolver3::computeAdvection(timeIntervalInSeconds);
    advect(timeIntervalInSeconds);
}

void GridSmokeSolver3::computeDiffusion(double timeIntervalInSeconds) {
    if (diffusionSolver() != nullptr) {
        if (_smokeDiffusionCoefficient > vox::kEpsilonD) {
            auto den = smokeDensity();
            auto den0 = std::dynamic_pointer_cast<CellCenteredScalarGrid3>(
                                                                           den->clone());
            
            diffusionSolver()->solve(*den0, _smokeDiffusionCoefficient,
                                     timeIntervalInSeconds, den.get(),
                                     *colliderSdf());
            extrapolateIntoCollider(den.get());
        }
        
        if (_temperatureDiffusionCoefficient > vox::kEpsilonD) {
            auto temp = smokeDensity();
            auto temp0 = std::dynamic_pointer_cast<CellCenteredScalarGrid3>(
                                                                            temp->clone());
            
            diffusionSolver()->solve(*temp0, _temperatureDiffusionCoefficient,
                                     timeIntervalInSeconds, temp.get(),
                                     *colliderSdf());
            extrapolateIntoCollider(temp.get());
        }
    }
    
    auto den = smokeDensity();
    den->parallelForEachDataPointIndex([&](const openvdb::Coord& coord) {
        (*den).getGrid()->tree().setValue(coord,
                                          (*den).getGrid()->tree().getValue(coord)
                                          *( 1.0 - _smokeDecayFactor));
    });
    auto temp = temperature();
    temp->parallelForEachDataPointIndex([&](const openvdb::Coord& coord) {
        (*temp).getGrid()->tree().setValue(coord,
                                           (*temp).getGrid()->tree().getValue(coord)
                                           *( 1.0 - _temperatureDecayFactor));
    });
}

void GridSmokeSolver3::computeBuoyancyForce(double timeIntervalInSeconds) {
    auto grids = vdbSystemData();
    auto vel = grids->velocity();
    
    vox::Vector3D up(0, 1, 0);
    if (gravity().lengthSquared() > vox::kEpsilonD) {
        up = -gravity().normalized();
    }
    
    if (std::abs(_buoyancySmokeDensityFactor) > vox::kEpsilonD ||
        std::abs(_buoyancyTemperatureFactor) > vox::kEpsilonD) {
        auto den = smokeDensity();
        auto temp = temperature();
        
        double tAmb = 0.0;
        temp->forEachCellIndex(
                               [&](uint i, uint j, uint k) {
            tAmb += (*temp)(openvdb::Coord(i, j, k));
        });
        tAmb /= static_cast<double>(
                                    temp->resolution().x
                                    * temp->resolution().y
                                    * temp->resolution().z);
        
        auto uPos = vel->uPosition();
        auto vPos = vel->vPosition();
        auto wPos = vel->wPosition();
        
        if (std::abs(up.x) > vox::kEpsilonD) {
            vel->parallelForEachUIndex([&](const openvdb::Coord& coord) {
                vox::Vector3D pt = uPos(coord);
                double fBuoy =
                _buoyancySmokeDensityFactor * den->sample(pt) +
                _buoyancyTemperatureFactor * (temp->sample(pt) - tAmb);
                vel->getUGrid()->tree().setValue(coord,
                                                 vel->getUGrid()->tree().getValue(coord)
                                                 + timeIntervalInSeconds * fBuoy * up.x);
            });
        }
        
        if (std::abs(up.y) > vox::kEpsilonD) {
            vel->parallelForEachVIndex([&](const openvdb::Coord& coord) {
                vox::Vector3D pt = vPos(coord);
                double fBuoy =
                _buoyancySmokeDensityFactor * den->sample(pt) +
                _buoyancyTemperatureFactor * (temp->sample(pt) - tAmb);
                vel->getVGrid()->tree().setValue(coord,
                                                 vel->getVGrid()->tree().getValue(coord)
                                                 + timeIntervalInSeconds * fBuoy * up.y);
            });
        }
        
        if (std::abs(up.z) > vox::kEpsilonD) {
            vel->parallelForEachWIndex([&](const openvdb::Coord& coord) {
                vox::Vector3D pt = wPos(coord);
                double fBuoy =
                _buoyancySmokeDensityFactor * den->sample(pt) +
                _buoyancyTemperatureFactor * (temp->sample(pt) - tAmb);
                vel->getWGrid()->tree().setValue(coord,
                                                 vel->getWGrid()->tree().getValue(coord)
                                                 + timeIntervalInSeconds * fBuoy * up.z);
            });
        }
        
        applyBoundaryCondition();
    }
}

GridSmokeSolver3::Builder GridSmokeSolver3::builder() { return Builder(); }

GridSmokeSolver3 GridSmokeSolver3::Builder::build() const {
    return GridSmokeSolver3(_resolution, getGridSpacing(), _gridOrigin);
}

GridSmokeSolver3Ptr GridSmokeSolver3::Builder::makeShared() const {
    return std::shared_ptr<GridSmokeSolver3>(
                                             new GridSmokeSolver3(_resolution,
                                                                  getGridSpacing(),
                                                                  _gridOrigin),
                                             [](GridSmokeSolver3* obj) { delete obj; });
}
