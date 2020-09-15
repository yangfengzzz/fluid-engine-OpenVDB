// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/array_utils.h"
#include "../src.common/fdm_utils.h"
#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_eno_level_set_solver3.hpp"
#include "../src.vdb/vdb_fmm_level_set_solver3.hpp"
#include "../src.vdb/vdb_upwind_level_set_solver3.hpp"

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(UpwindLevelSetSolver3, Reinitialize) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);

    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);

    UpwindLevelSetSolver3 solver;
    solver.reinitialize(sdf, 5.0, &temp);

    for (uint k = 0; k < 50; ++k) {
        for (uint j = 0; j < 30; ++j) {
            for (uint i = 0; i < 40; ++i) {
                EXPECT_NEAR(sdf(openvdb::Coord(i, j, k)),
                            temp(openvdb::Coord(i, j, k)), 0.7)
                << i << ", " << j << ", " << k;
            }
        }
    }
}

TEST(UpwindLevelSetSolver3, Extrapolate) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);
    CellCenteredScalarGrid3 field(40, 30, 50);

    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);
    field.fill(5.0, vox::ExecutionPolicy::kSerial);

    UpwindLevelSetSolver3 solver;
    solver.extrapolate(field, sdf, 5.0, &temp);

    for (uint k = 0; k < 50; ++k) {
        for (uint j = 0; j < 30; ++j) {
            for (uint i = 0; i < 40; ++i) {
                EXPECT_DOUBLE_EQ(5.0, temp(openvdb::Coord(i, j, k) ))
                << i << ", " << j << ", " << k;
            }
        }
    }
}


//--------------------------------------------------------
TEST(EnoLevelSetSolver3, Reinitialize) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);
    
    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);

    EnoLevelSetSolver3 solver;
    solver.reinitialize(sdf, 5.0, &temp);

    for (uint k = 0; k < 50; ++k) {
        for (uint j = 0; j < 30; ++j) {
            for (uint i = 0; i < 40; ++i) {
                EXPECT_NEAR(sdf(openvdb::Coord(i, j, k)),
                            temp(openvdb::Coord(i, j, k)), 0.5)
                << i << ", " << j << ", " << k;
            }
        }
    }
}

TEST(EnoLevelSetSolver3, Extrapolate) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);
    CellCenteredScalarGrid3 field(40, 30, 50);

    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);
    field.fill(5.0, vox::ExecutionPolicy::kSerial);

    EnoLevelSetSolver3 solver;
    solver.extrapolate(field, sdf, 5.0, &temp);

    for (uint k = 0; k < 50; ++k) {
        for (uint j = 0; j < 30; ++j) {
            for (uint i = 0; i < 40; ++i) {
                EXPECT_DOUBLE_EQ(5.0, temp(openvdb::Coord(i, j, k)))
                << i << ", " << j << ", " << k;
            }
        }
    }
}

//------------------------------------------------------------
TEST(FmmLevelSetSolver3, Reinitialize) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);
    
    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);
    
    FmmLevelSetSolver3 solver;
    solver.reinitialize(sdf, 5.0, &temp);
    
    for (uint k = 0; k < 50; ++k) {
        for (uint j = 0; j < 30; ++j) {
            for (uint i = 0; i < 40; ++i) {
                EXPECT_NEAR(sdf(openvdb::Coord(i, j, k)),
                            temp(openvdb::Coord(i, j, k)), 0.9)
                << i << ", " << j << ", " << k;
            }
        }
    }
}

TEST(FmmLevelSetSolver3, Extrapolate) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);
    CellCenteredScalarGrid3 field(40, 30, 50);

    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);
    field.fill(5.0, vox::ExecutionPolicy::kSerial);

    FmmLevelSetSolver3 solver;
    solver.extrapolate(field, sdf, 5.0, &temp);

    for (uint k = 0; k < 50; ++k) {
        for (uint j = 0; j < 30; ++j) {
            for (uint i = 0; i < 40; ++i) {
                EXPECT_DOUBLE_EQ(5.0, temp(openvdb::Coord(i, j, k)) )
                << i << ", " << j << ", " << k;
            }
        }
    }
}
