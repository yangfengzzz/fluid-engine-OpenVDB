//
//  fdm_utils_tests.cpp
//  vdb_tests
//
//  Created by Feng Yang on 2020/2/20.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_cell_centered_vector_grid3.h"
#include "../src.vdb/vdb_fdm_utils.h"

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(FdmUtils, ScalarToGradient3) {
    CellCenteredScalarGrid3 grid(10, 10, 10, 2.0, 3.0, 0.5, -1.0, 4.0, 2.0);
    grid.fill([&](const vox::Vector3D& x) {
        return -5.0 * x.x + 4.0 * x.y + 2.0 * x.z;
    }, vox::ExecutionPolicy::kSerial);
    
    vox::Vector3D grad = gradient3(
                                   grid.getGrid(),
                                   grid.dataSize(),
                                   grid.gridSpacing(),
                                   5, 3, 4);
    EXPECT_DOUBLE_EQ(-5.0, grad.x);
    EXPECT_DOUBLE_EQ(4.0, grad.y);
    EXPECT_DOUBLE_EQ(2.0, grad.z);
}

TEST(FdmUtils, VectorToGradient3) {
    CellCenteredVectorGrid3 grid(10, 10, 10, 2.0, 3.0, 0.5, -1.0, 4.0, 2.0);
    grid.fill([&](const vox::Vector3D& x) {
        return vox::Vector3D(
                             -5.0 * x.x + 4.0 * x.y + 2.0 * x.z,
                             2.0 * x.x - 7.0 * x.y,
                             x.y + 3.0 * x.z);
    }, vox::ExecutionPolicy::kSerial);
    
    auto grad = gradient3(
                          grid.getGrid(),
                          grid.dataSize(),
                          grid.gridSpacing(),
                          5, 3, 4);
    EXPECT_DOUBLE_EQ(-5.0, grad[0].x);
    EXPECT_DOUBLE_EQ(4.0, grad[0].y);
    EXPECT_DOUBLE_EQ(2.0, grad[0].z);
    EXPECT_DOUBLE_EQ(2.0, grad[1].x);
    EXPECT_DOUBLE_EQ(-7.0, grad[1].y);
    EXPECT_DOUBLE_EQ(0.0, grad[1].z);
    EXPECT_DOUBLE_EQ(0.0, grad[2].x);
    EXPECT_DOUBLE_EQ(1.0, grad[2].y);
    EXPECT_DOUBLE_EQ(3.0, grad[2].z);
}

TEST(FdmUtils, ScalarToLaplacian3) {
    CellCenteredScalarGrid3 grid(10, 10, 10, 2.0, 3.0, 0.5, -1.0, 4.0, 2.0);
    grid.fill([&](const vox::Vector3D& x) {
        return -5.0 * x.x * x.x + 4.0 * x.y * x.y - 3.0 * x.z * x.z;
    }, vox::ExecutionPolicy::kSerial);
    
    double lapl = laplacian3(
                             grid.getGrid(),
                             grid.dataSize(),
                             grid.gridSpacing(),
                             5, 3, 4);
    EXPECT_DOUBLE_EQ(-8.0, lapl);
}

TEST(FdmUtils, VectorToLaplacian3) {
    CellCenteredVectorGrid3 grid(10, 10, 10, 2.0, 3.0, 0.5, -1.0, 4.0, 2.0);
    grid.fill([&](const vox::Vector3D& x) {
        return vox::Vector3D(
                             -5.0 * x.x * x.x + 4.0 * x.y * x.y + 2.0 * x.z * x.z,
                             2.0 * x.x * x.x - 7.0 * x.y * x.y,
                             x.y * x.y + 3.0 * x.z * x.z);
    }, vox::ExecutionPolicy::kSerial);
    
    auto lapl = laplacian3(
                           grid.getGrid(),
                           grid.dataSize(),
                           grid.gridSpacing(),
                           5, 3, 4);
    EXPECT_DOUBLE_EQ(2.0, lapl.x);
    EXPECT_DOUBLE_EQ(-10.0, lapl.y);
    EXPECT_DOUBLE_EQ(8.0, lapl.z);
}
