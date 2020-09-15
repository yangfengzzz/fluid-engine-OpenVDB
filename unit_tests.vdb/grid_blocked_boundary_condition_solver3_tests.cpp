// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.vdb/vdb_blocked_boundary_condition_solver3.hpp"
#include "../src.common/plane3.h"
#include "../src.common/rigid_body_collider3.h"

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(GridBlockedBoundaryConditionSolver3, ClosedDomain) {
    BlockedBoundaryConditionSolver3 bndSolver;
    vox::Size3 gridSize(10, 10, 10);
    vox::Vector3D gridSpacing(1.0, 1.0, 1.0);
    vox::Vector3D gridOrigin(-5.0, -5.0, -5.0);
    
    bndSolver.updateCollider(nullptr, gridSize, gridSpacing, gridOrigin);
    
    FaceCenteredGrid3 velocity(gridSize, gridSpacing, gridOrigin);
    velocity.fill(vox::Vector3D(1.0, 1.0, 1.0), vox::ExecutionPolicy::kSerial);
    
    bndSolver.constrainVelocity(&velocity);
    
    velocity.forEachUIndex([&](const openvdb::Coord& coord) {
        if (coord.x() == 0 || coord.x() == gridSize.x) {
            EXPECT_DOUBLE_EQ(0.0, velocity.u(coord));
        } else {
            EXPECT_DOUBLE_EQ(1.0, velocity.u(coord));
        }
    });
    
    velocity.forEachVIndex([&](const openvdb::Coord& coord) {
        if (coord.y() == 0 || coord.y() == gridSize.y) {
            EXPECT_DOUBLE_EQ(0.0, velocity.v(coord));
        } else {
            EXPECT_DOUBLE_EQ(1.0, velocity.v(coord));
        }
    });
    
    velocity.forEachWIndex([&](const openvdb::Coord& coord) {
        if (coord.z() == 0 || coord.z() == gridSize.z) {
            EXPECT_DOUBLE_EQ(0.0, velocity.w(coord));
        } else {
            EXPECT_DOUBLE_EQ(1.0, velocity.w(coord));
        }
    });
}

TEST(GridBlockedBoundaryConditionSolver3, OpenDomain) {
    BlockedBoundaryConditionSolver3 bndSolver;
    vox::Size3 gridSize(10, 10, 10);
    vox::Vector3D gridSpacing(1.0, 1.0, 1.0);
    vox::Vector3D gridOrigin(-5.0, -5.0, -5.0);
    
    // Partially open domain
    bndSolver.setClosedDomainBoundaryFlag(
                                          vox::kDirectionLeft |
                                          vox::kDirectionUp |
                                          vox::kDirectionFront);
    bndSolver.updateCollider(nullptr, gridSize, gridSpacing, gridOrigin);
    
    FaceCenteredGrid3 velocity(gridSize, gridSpacing, gridOrigin);
    velocity.fill(vox::Vector3D(1.0, 1.0, 1.0), vox::ExecutionPolicy::kSerial);
    
    bndSolver.constrainVelocity(&velocity);
    
    velocity.forEachUIndex([&](const openvdb::Coord& coord) {
        if (coord.x() == 0) {
            EXPECT_DOUBLE_EQ(0.0, velocity.u(coord));
        } else {
            EXPECT_DOUBLE_EQ(1.0, velocity.u(coord));
        }
    });
    
    velocity.forEachVIndex([&](const openvdb::Coord& coord) {
        if (coord.y() == gridSize.y) {
            EXPECT_DOUBLE_EQ(0.0, velocity.v(coord));
        } else {
            EXPECT_DOUBLE_EQ(1.0, velocity.v(coord));
        }
    });
    
    velocity.forEachWIndex([&](const openvdb::Coord& coord) {
        if (coord.z() == gridSize.z) {
            EXPECT_DOUBLE_EQ(0.0, velocity.w(coord));
        } else {
            EXPECT_DOUBLE_EQ(1.0, velocity.w(coord));
        }
    });
}
