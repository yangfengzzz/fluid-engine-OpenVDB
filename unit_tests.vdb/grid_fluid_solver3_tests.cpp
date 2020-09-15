// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.vdb/vdb_fluid_solver3.hpp"

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(GridFluidSolver3, Constructor) {
    FluidSolver3 solver;
    
    // Check if the sub-step solvers are present
    EXPECT_TRUE(solver.diffusionSolver() != nullptr);
    EXPECT_TRUE(solver.pressureSolver() != nullptr);
    
    // Check default parameters
    EXPECT_GE(solver.viscosityCoefficient(), 0.0);
    EXPECT_GT(solver.maxCfl(), 0.0);
    EXPECT_EQ(vox::kDirectionAll, solver.closedDomainBoundaryFlag());
    
    // Check grid system data
    EXPECT_TRUE(solver.vdbSystemData() != nullptr);
    EXPECT_EQ(1u, solver.vdbSystemData()->resolution().x);
    EXPECT_EQ(1u, solver.vdbSystemData()->resolution().y);
    EXPECT_EQ(1u, solver.vdbSystemData()->resolution().z);
    EXPECT_EQ(solver.vdbSystemData()->velocity(), solver.velocity());
    
    // Collider should be null
    EXPECT_TRUE(solver.collider() == nullptr);
}

TEST(GridFluidSolver3, ResizeGridSystemData) {
    FluidSolver3 solver;
    
    solver.resizeGrid(
                      vox::Size3(1, 2, 3),
                      vox::Vector3D(4.0, 5.0, 6.0),
                      vox::Vector3D(7.0, 8.0, 9.0));
    
    EXPECT_EQ(1u, solver.resolution().x);
    EXPECT_EQ(2u, solver.resolution().y);
    EXPECT_EQ(3u, solver.resolution().z);
    EXPECT_EQ(4.0, solver.gridSpacing().x);
    EXPECT_EQ(5.0, solver.gridSpacing().y);
    EXPECT_EQ(6.0, solver.gridSpacing().z);
    EXPECT_EQ(7.0, solver.gridOrigin().x);
    EXPECT_EQ(8.0, solver.gridOrigin().y);
    EXPECT_EQ(9.0, solver.gridOrigin().z);
}

TEST(GridFluidSolver3, MinimumResolution) {
    FluidSolver3 solver;
    
    solver.resizeGrid(vox::Size3(1, 1, 1), vox::Vector3D(1.0, 1.0, 1.0), vox::Vector3D());
    solver.velocity()->fill(vox::Vector3D(), vox::ExecutionPolicy::kSerial);
    
    vox::Frame frame(0, 1.0 / 60.0);
    frame.timeIntervalInSeconds = 0.01;
    solver.update(frame);
}

TEST(GridFluidSolver3, GravityOnly) {
    FluidSolver3 solver;
    solver.setGravity(vox::Vector3D(0, -10, 0.0));
    solver.setDiffusionSolver(nullptr);
    solver.setPressureSolver(nullptr);
    
    solver.resizeGrid(
                      vox::Size3(3, 3, 3),
                      vox::Vector3D(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
                      vox::Vector3D());
    solver.velocity()->fill(vox::Vector3D(), vox::ExecutionPolicy::kSerial);
    
    vox::Frame frame(0, 0.01);
    solver.update(frame);
    
    solver.velocity()->forEachUIndex([&](const openvdb::Coord& coord) {
        EXPECT_NEAR(0.0, solver.velocity()->u(coord), 1e-8);
    });
    
    solver.velocity()->forEachVIndex([&](const openvdb::Coord& coord) {
        if (coord.y() == 0 || coord.y() == 3) {
            EXPECT_NEAR(0.0, solver.velocity()->v(coord), 1e-8);
        } else {
            EXPECT_NEAR(-0.1, solver.velocity()->v(coord), 1e-8);
        }
    });
    
    solver.velocity()->forEachWIndex([&](const openvdb::Coord& coord) {
        EXPECT_NEAR(0.0, solver.velocity()->w(coord), 1e-8);
    });
}
