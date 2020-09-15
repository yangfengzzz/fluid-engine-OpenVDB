// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_backward_euler_diffusion_solver3.hpp"
#include "../src.common/array3.h"

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(GridBackwardEulerDiffusionSolver3, Solve) {
    CellCenteredScalarGrid3 src(3, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    CellCenteredScalarGrid3 dst(3, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    
    src.dataAccessor().setValue(openvdb::Coord(1, 1, 1), 1.0);
    
    BackwardEulerDiffusionSolver3 diffusionSolver;
    diffusionSolver.solve(src, 1.0 / 12.0, 1.0, &dst);
    
    vox::Array3<double> solution = {
        {
            {0.001058, 0.005291, 0.001058},
            {0.005291, 0.041270, 0.005291},
            {0.001058, 0.005291, 0.001058}
        },
        {
            {0.005291, 0.041270, 0.005291},
            {0.041270, 0.680423, 0.041270},
            {0.005291, 0.041270, 0.005291}
        },
        {
            {0.001058, 0.005291, 0.001058},
            {0.005291, 0.041270, 0.005291},
            {0.001058, 0.005291, 0.001058}
        }
    };
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                EXPECT_NEAR(solution(i, j, k),
                            dst(openvdb::Coord(i, j, k)), 1e-6);
            }
        }
    }
}
