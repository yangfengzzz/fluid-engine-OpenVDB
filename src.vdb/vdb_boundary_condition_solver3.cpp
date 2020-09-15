//
//  vdb_boundary_condition_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_boundary_condition_solver3.hpp"

using namespace vdb;

BoundaryConditionSolver3::BoundaryConditionSolver3() {
}

BoundaryConditionSolver3::~BoundaryConditionSolver3() {
}

const vox::Collider3Ptr& BoundaryConditionSolver3::collider() const {
    return _collider;
}

void BoundaryConditionSolver3::updateCollider(
                                              const vox::Collider3Ptr& newCollider,
                                              const vox::Size3& gridSize,
                                              const vox::Vector3D& gridSpacing,
                                              const vox::Vector3D& gridOrigin) {
    _collider = newCollider;
    _gridSize = gridSize;
    _gridSpacing = gridSpacing;
    _gridOrigin = gridOrigin;
    
    onColliderUpdated(gridSize, gridSpacing, gridOrigin);
}

int BoundaryConditionSolver3::closedDomainBoundaryFlag() const {
    return _closedDomainBoundaryFlag;
}

void BoundaryConditionSolver3::setClosedDomainBoundaryFlag(int flag) {
    _closedDomainBoundaryFlag = flag;
}

const vox::Size3& BoundaryConditionSolver3::gridSize() const {
    return _gridSize;
}

const vox::Vector3D& BoundaryConditionSolver3::gridSpacing() const {
    return _gridSpacing;
}

const vox::Vector3D& BoundaryConditionSolver3::gridOrigin() const {
    return _gridOrigin;
}
