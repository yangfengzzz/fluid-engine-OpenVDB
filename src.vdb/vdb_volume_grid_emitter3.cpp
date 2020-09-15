//
//  _grid_emitter3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/6.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"

#include "../src.common/level_set_utils.h"
#include "vdb_collocated_vector_grid3.h"
#include "vdb_face_centered_grid3.h"
#include "../src.common/surface_to_implicit3.h"
#include "vdb_volume_grid_emitter3.hpp"

#include <algorithm>

using namespace vdb;

GridEmitter3::GridEmitter3(const vox::ImplicitSurface3Ptr& sourceRegion,
                           bool isOneShot)
: _sourceRegion(sourceRegion), _isOneShot(isOneShot) {}

GridEmitter3::~GridEmitter3() {}

void GridEmitter3::addSignedDistanceTarget(
                                           const ScalarGrid3Ptr& scalarGridTarget) {
    auto mapper = [](double sdf, const vox::Vector3D&, double oldVal) {
        return std::min(oldVal, sdf);
    };
    addTarget(scalarGridTarget, mapper);
}

void GridEmitter3::addStepFunctionTarget(
                                         const ScalarGrid3Ptr& scalarGridTarget,
                                         double minValue, double maxValue) {
    double smoothingWidth = scalarGridTarget->gridSpacing().min();
    auto mapper = [minValue, maxValue, smoothingWidth, scalarGridTarget](
                                                                         double sdf,
                                                                         const vox::Vector3D&,
                                                                         double oldVal) {
        double step = 1.0 - vox::smearedHeavisideSdf(sdf / smoothingWidth);
        return std::max(oldVal, (maxValue - minValue) * step + minValue);
    };
    addTarget(scalarGridTarget, mapper);
}

void GridEmitter3::addTarget(const ScalarGrid3Ptr& scalarGridTarget,
                             const ScalarMapper& customMapper) {
    _customScalarTargets.emplace_back(scalarGridTarget, customMapper);
}

void GridEmitter3::addTarget(const VectorGrid3Ptr& vectorGridTarget,
                             const VectorMapper& customMapper) {
    _customVectorTargets.emplace_back(vectorGridTarget, customMapper);
}

const vox::ImplicitSurface3Ptr& GridEmitter3::sourceRegion() const {
    return _sourceRegion;
}

bool GridEmitter3::isOneShot() const { return _isOneShot; }

//---------------------------------------------------------------
void GridEmitter3::onUpdate(double currentTimeInSeconds,
                            double timeIntervalInSeconds) {
    UNUSED_VARIABLE(currentTimeInSeconds);
    UNUSED_VARIABLE(timeIntervalInSeconds);
    
    if (!isEnabled()) {
        return;
    }
    
    emit();
    
    if (_isOneShot) {
        setIsEnabled(false);
    }
    
    _hasEmitted = true;
}

void GridEmitter3::emit() {
    if (!_sourceRegion) {
        return;
    }
    
    _sourceRegion->updateQueryEngine();
    
    for (const auto& target : _customScalarTargets) {
        const auto& grid = std::get<0>(target);
        const auto& mapper = std::get<1>(target);
        
        auto pos = grid->dataPosition();
        grid->parallelForEachDataPointIndex([&](const openvdb::Coord& coord) {
            vox::Vector3D gx = pos(coord);
            double sdf = sourceRegion()->signedDistance(gx);
            grid->getGrid()->tree().setValue(coord, mapper(sdf, gx, grid->getGrid()->tree().getValue(coord)) );
        });
    }
    
    for (const auto& target : _customVectorTargets) {
        const auto& grid = std::get<0>(target);
        const auto& mapper = std::get<1>(target);
        
        CollocatedVectorGrid3Ptr collocated =
        std::dynamic_pointer_cast<CollocatedVectorGrid3>(grid);
        if (collocated != nullptr) {
            auto pos = collocated->dataPosition();
            collocated->parallelForEachDataPointIndex(
                                                      [&](const openvdb::Coord& coord) {
                vox::Vector3D gx = pos(coord);
                double sdf = sourceRegion()->signedDistance(gx);
                if (vox::isInsideSdf(sdf)) {
                    openvdb::Vec3d val = collocated->getGrid()->tree().getValue(coord);
                    vox::Vector3D m_val = mapper(sdf, gx,
                                                 vox::Vector3D(val.x(), val.y(), val.z()));
                    collocated->getGrid()->tree().setValue(coord,
                                                           openvdb::Vec3d(m_val.x, m_val.y, m_val.z) );
                }
            });
            continue;
        }
        
        FaceCenteredGrid3Ptr faceCentered =
        std::dynamic_pointer_cast<FaceCenteredGrid3>(grid);
        if (faceCentered != nullptr) {
            auto uPos = faceCentered->uPosition();
            auto vPos = faceCentered->vPosition();
            auto wPos = faceCentered->wPosition();
            
            faceCentered->parallelForEachUIndex(
                                                [&](const openvdb::Coord& coord) {
                vox::Vector3D gx = uPos(coord);
                double sdf = sourceRegion()->signedDistance(gx);
                vox::Vector3D oldVal = faceCentered->sample(gx);
                vox::Vector3D newVal = mapper(sdf, gx, oldVal);
                faceCentered->getUGrid()->tree().setValue(coord, newVal.x);
            });
            faceCentered->parallelForEachVIndex(
                                                [&](const openvdb::Coord& coord) {
                vox::Vector3D gx = vPos(coord);
                double sdf = sourceRegion()->signedDistance(gx);
                vox::Vector3D oldVal = faceCentered->sample(gx);
                vox::Vector3D newVal = mapper(sdf, gx, oldVal);
                faceCentered->getVGrid()->tree().setValue(coord, newVal.y);
            });
            faceCentered->parallelForEachWIndex(
                                                [&](const openvdb::Coord& coord) {
                vox::Vector3D gx = wPos(coord);
                double sdf = sourceRegion()->signedDistance(gx);
                vox::Vector3D oldVal = faceCentered->sample(gx);
                vox::Vector3D newVal = mapper(sdf, gx, oldVal);
                faceCentered->getWGrid()->tree().setValue(coord, newVal.z);
            });
            continue;
        }
    }
}

//---------------------------------------------------------------
GridEmitter3::Builder GridEmitter3::builder() { return Builder(); }

GridEmitter3::Builder& GridEmitter3::Builder::withSourceRegion(
                                                               const vox::Surface3Ptr& sourceRegion) {
    auto implicit = std::dynamic_pointer_cast<vox::ImplicitSurface3>(sourceRegion);
    if (implicit != nullptr) {
        _sourceRegion = implicit;
    } else {
        _sourceRegion = std::make_shared<vox::SurfaceToImplicit3>(sourceRegion);
    }
    return *this;
}

GridEmitter3::Builder& GridEmitter3::Builder::withIsOneShot(
                                                            bool isOneShot) {
    _isOneShot = isOneShot;
    return *this;
}

GridEmitter3 GridEmitter3::Builder::build() const {
    return GridEmitter3(_sourceRegion, _isOneShot);
}

GridEmitter3Ptr GridEmitter3::Builder::makeShared() const {
    return std::shared_ptr<GridEmitter3>(
                                         new GridEmitter3(_sourceRegion, _isOneShot),
                                         [](GridEmitter3* obj) { delete obj; });
}
