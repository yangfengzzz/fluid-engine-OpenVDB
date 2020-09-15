// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_forward_euler_diffusion_solver3.hpp"

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(GridForwardEulerDiffusionSolver3, Solve) {
    CellCenteredScalarGrid3 src(3, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    CellCenteredScalarGrid3 dst(3, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

    src.getGrid()->getAccessor().setValue(openvdb::Coord(1, 1, 1), 1.0);
    
    ForwardEulerDiffusionSolver3 diffusionSolver;
    diffusionSolver.solve(src, 1.0 / 12.0, 1.0, &dst);
    
    EXPECT_DOUBLE_EQ(1.0/12.0, dst.getGrid()->getAccessor().getValue(openvdb::Coord(1, 1, 0)));
    EXPECT_DOUBLE_EQ(1.0/12.0, dst.getGrid()->getAccessor().getValue(openvdb::Coord(0, 1, 1)));
    EXPECT_DOUBLE_EQ(1.0/12.0, dst.getGrid()->getAccessor().getValue(openvdb::Coord(1, 0, 1)));
    EXPECT_DOUBLE_EQ(1.0/12.0, dst.getGrid()->getAccessor().getValue(openvdb::Coord(2, 1, 1)));
    EXPECT_DOUBLE_EQ(1.0/12.0, dst.getGrid()->getAccessor().getValue(openvdb::Coord(1, 2, 1)));
    EXPECT_DOUBLE_EQ(1.0/12.0, dst.getGrid()->getAccessor().getValue(openvdb::Coord(1, 1, 2)));
    EXPECT_DOUBLE_EQ(1.0/2.0,  dst.getGrid()->getAccessor().getValue(openvdb::Coord(1, 1, 1)));
}
