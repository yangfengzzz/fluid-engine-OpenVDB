//
//  multiRes_tests.cpp
//  vdb_tests
//
//  Created by Feng Yang on 2020/2/10.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_face_centered_grid3.h"
#include "../src.vdb/vdb_helper.h"
#include <openvdb/tools/MultiResGrid.h>

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(MultiResolution, MultiResGridConstruction) {
    const double background = -1.0;
    const size_t levels = 4;
    
    openvdb::tools::MultiResGrid<openvdb::DoubleTree>::Ptr
    weights(new openvdb::tools::MultiResGrid<openvdb::DoubleTree>(levels, background));
    EXPECT_EQ(levels  , weights->numLevels());
    EXPECT_EQ(size_t(0), weights->finestLevel());
    EXPECT_EQ(levels-1, weights->coarsestLevel());
}

TEST(MultiResolution, MultiResGridConstruction2) {
    openvdb::DoubleGrid::Ptr finest = openvdb::DoubleGrid::create(0.0);
    openvdb::math::ScaleTranslateMap transMap(openvdb::Vec3d(0.001, 0.001, 0.001),
                                              openvdb::Vec3d(1, 2, 3));
    openvdb::math::Transform trans(transMap.copy());
    finest->setTransform(trans.copy());
    for (int i = 0; i < 200; i++) {
        for (int j = 0; j < 200; j++) {
            for (int k = 0; k < 200; k++) {
                finest->tree().setValueOn(openvdb::Coord(i, j, k), 100);
            }
        }
    }
    const size_t levels = 4;
    
    openvdb::tools::MultiResGrid<openvdb::DoubleTree>::Ptr
    weights(new openvdb::tools::MultiResGrid<openvdb::DoubleTree>(levels, finest,
                                                                  true));//must use injection
    EXPECT_EQ(levels  , weights->numLevels());
    EXPECT_EQ(size_t(0), weights->finestLevel());
    EXPECT_EQ(levels-1, weights->coarsestLevel());
    
    EXPECT_EQ(weights->grid(0)->activeVoxelCount(), 200*200*200);
    EXPECT_EQ(weights->grid(1)->activeVoxelCount(), 100*100*100);
    EXPECT_EQ(weights->grid(2)->activeVoxelCount(), 50*50*50);
    EXPECT_EQ(weights->grid(3)->activeVoxelCount(), 25*25*25);
        
    EXPECT_NEAR(weights->grid(0)->voxelSize().x(), 0.001, 1.0e-15);
    EXPECT_NEAR(weights->grid(0)->voxelSize().y(), 0.001, 1.0e-15);
    EXPECT_NEAR(weights->grid(0)->voxelSize().z(), 0.001, 1.0e-15);
    EXPECT_NEAR(weights->grid(1)->voxelSize().x(), 0.002, 1.0e-15);
    EXPECT_NEAR(weights->grid(1)->voxelSize().y(), 0.002, 1.0e-15);
    EXPECT_NEAR(weights->grid(1)->voxelSize().z(), 0.002, 1.0e-15);
    EXPECT_NEAR(weights->grid(2)->voxelSize().x(), 0.004, 1.0e-15);
    EXPECT_NEAR(weights->grid(2)->voxelSize().y(), 0.004, 1.0e-15);
    EXPECT_NEAR(weights->grid(2)->voxelSize().z(), 0.004, 1.0e-15);
    EXPECT_NEAR(weights->grid(3)->voxelSize().x(), 0.008, 1.0e-15);
    EXPECT_NEAR(weights->grid(3)->voxelSize().y(), 0.008, 1.0e-15);
    EXPECT_NEAR(weights->grid(3)->voxelSize().z(), 0.008, 1.0e-15);
    
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 25; j++) {
            for (int k = 0; k < 25; k++) {
                EXPECT_EQ(weights->sampleValue<1>(openvdb::Coord(i, j, k), 3), 100);
            }
        }
    }
    
    EXPECT_EQ(weights->grid(0)->indexToWorld(openvdb::Coord(0, 0, 0)), openvdb::Vec3d(1, 2, 3));
    EXPECT_EQ(weights->grid(1)->indexToWorld(openvdb::Coord(0, 0, 0)), openvdb::Vec3d(1, 2, 3));
    EXPECT_EQ(weights->grid(2)->indexToWorld(openvdb::Coord(0, 0, 0)), openvdb::Vec3d(1, 2, 3));
    EXPECT_EQ(weights->grid(3)->indexToWorld(openvdb::Coord(0, 0, 0)), openvdb::Vec3d(1, 2, 3));
}

TEST(MultiResolution, CoordBoxShift) {
    openvdb::CoordBBox box(openvdb::Coord(0, 0, 0), openvdb::Coord(200, 200, 200));
    box = box >> (size_t)1;
    EXPECT_EQ(box.min(), openvdb::Coord(0, 0, 0));
    EXPECT_EQ(box.max(), openvdb::Coord(100, 100, 100));
    box = box << (size_t)1;
    EXPECT_EQ(box.min(), openvdb::Coord(0, 0, 0));
    EXPECT_EQ(box.max(), openvdb::Coord(200, 200, 200));
}
