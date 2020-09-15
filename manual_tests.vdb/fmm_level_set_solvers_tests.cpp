// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "manual_tests.h"

#include "../src.common/array_utils.h"
#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_fmm_level_set_solver3.hpp"

using namespace vdb;

JET_TESTS(FmmLevelSetSolver3);

JET_BEGIN_TEST_F(FmmLevelSetSolver3, ReinitializeSmall) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);
    
    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);
    
    FmmLevelSetSolver3 solver;
    solver.reinitialize(sdf, 5.0, &temp);
    
    vox::Array2<double> sdf2(40, 30);
    vox::Array2<double> temp2(40, 30);
    for (uint j = 0; j < 30; ++j) {
        for (uint i = 0; i < 40; ++i) {
            sdf2(i, j) = sdf.getGrid()->tree().getValue(openvdb::Coord(i, j, 10));
            temp2(i, j) = temp.getGrid()->tree().getValue(openvdb::Coord(i, j, 10));
        }
    }
    
    saveData(sdf2.constAccessor(), "sdf_#grid2,iso.npy");
    saveData(temp2.constAccessor(), "temp_#grid2,iso.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(FmmLevelSetSolver3, ExtrapolateSmall) {
    CellCenteredScalarGrid3 sdf(40, 30, 50), temp(40, 30, 50);
    CellCenteredScalarGrid3 field(40, 30, 50);
    
    sdf.fill([](const vox::Vector3D& x) {
        return (x - vox::Vector3D(20, 20, 20)).length() - 8.0;
    }, vox::ExecutionPolicy::kSerial);
    field.fill([&](const vox::Vector3D& x) {
        if ((x - vox::Vector3D(20, 20, 20)).length() <= 8.0) {
            return 0.5 * 0.25 * std::sin(x.x) * std::sin(x.y) * std::sin(x.z);
        } else {
            return 0.0;
        }
    }, vox::ExecutionPolicy::kSerial);
    
    FmmLevelSetSolver3 solver;
    solver.extrapolate(field, sdf, 5.0, &temp);
    
    vox::Array2<double> field2(40, 30);
    vox::Array2<double> temp2(40, 30);
    for (uint j = 0; j < 30; ++j) {
        for (uint i = 0; i < 40; ++i) {
            field2(i, j) = field.getGrid()->tree().getValue(openvdb::Coord(i, j, 12));
            temp2(i, j) = temp.getGrid()->tree().getValue(openvdb::Coord(i, j, 12));
        }
    }
    
    saveData(field2.constAccessor(), "field_#grid2.npy");
    saveData(temp2.constAccessor(), "temp_#grid2.npy");
}
JET_END_TEST_F
