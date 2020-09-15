//
//  vdb_blocked_boundary_condition_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "../src.common/physics_helpers.h"
#include "../src.common/array_utils.h"
#include "vdb_blocked_boundary_condition_solver3.hpp"
#include "../src.common/level_set_utils.h"
#include "../src.common/surface_to_implicit3.h"
#include <algorithm>

using namespace vdb;

static const char kFluid = 1;
static const char kCollider = 0;

BlockedBoundaryConditionSolver3::BlockedBoundaryConditionSolver3() {
}

const vox::Array3<char>& BlockedBoundaryConditionSolver3::marker() const {
    return _marker;
}

void BlockedBoundaryConditionSolver3::onColliderUpdated(const vox::Size3& gridSize,
                                                           const vox::Vector3D& gridSpacing,
                                                           const vox::Vector3D& gridOrigin) {
    FractionalBoundaryConditionSolver3::onColliderUpdated(gridSize,
                                                             gridSpacing,
                                                             gridOrigin);
    
    const auto sdf
    = std::dynamic_pointer_cast<CellCenteredScalarGrid3>(colliderSdf());
    
    _marker.resize(gridSize);
    _marker.parallelForEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf((*sdf)(openvdb::Coord(i, j, k)))) {
            _marker(i, j, k) = kCollider;
        } else {
            _marker(i, j, k) = kFluid;
        }
    });
}

void BlockedBoundaryConditionSolver3::constrainVelocity(FaceCenteredGrid3* velocity,
                                                           unsigned int extrapolationDepth) {
    FractionalBoundaryConditionSolver3::constrainVelocity(
                                                             velocity,
                                                             extrapolationDepth);
    // No-flux: project the velocity at the marker interface
    vox::Size3 size = velocity->resolution();
    auto uPos = velocity->uPosition();
    auto vPos = velocity->vPosition();
    auto wPos = velocity->wPosition();
    
    const auto sdf
    = std::dynamic_pointer_cast<CellCenteredScalarGrid3>(colliderSdf());
    
    _marker.forEachIndex([&](uint i, uint j, uint k) {
        if (_marker(i, j, k) == kCollider) {
            openvdb::Coord coord(i, j, k);
            if (i > 0 && _marker(i - 1, j, k) == kFluid) {
                vox::Vector3D colliderVel = collider()->velocityAt(uPos(coord));
                velocity->getUGrid()->tree().setValueOnly(coord, colliderVel.x);
            }
            
            openvdb::Coord neigh = coord + openvdb::Coord(1, 0, 0);
            if (i < size.x - 1 && _marker(i + 1, j, k) == kFluid) {
                vox::Vector3D colliderVel = collider()->velocityAt(uPos(neigh));
                velocity->getUGrid()->tree().setValueOnly(neigh, colliderVel.x);
            }
            
            if (j > 0 && _marker(i, j - 1, k) == kFluid) {
                vox::Vector3D colliderVel = collider()->velocityAt(vPos(coord));
                velocity->getVGrid()->tree().setValueOnly(coord, colliderVel.y);
            }
            
            neigh = coord + openvdb::Coord(0, 1, 0);
            if (j < size.y - 1 && _marker(i, j + 1, k) == kFluid) {
                vox::Vector3D colliderVel = collider()->velocityAt(vPos(neigh));
                velocity->getVGrid()->tree().setValueOnly(neigh, colliderVel.y);
            }
            
            if (k > 0 && _marker(i, j, k - 1) == kFluid) {
                vox::Vector3D colliderVel = collider()->velocityAt(wPos(coord));
                velocity->getWGrid()->tree().setValueOnly(coord, colliderVel.z);
            }
            
            neigh = coord + openvdb::Coord(0, 0, 1);
            if (k < size.z - 1 && _marker(i, j, k + 1) == kFluid) {
                vox::Vector3D colliderVel = collider()->velocityAt(wPos(neigh));
                velocity->getWGrid()->tree().setValueOnly(neigh, colliderVel.z);
            }
        }
    });
}
