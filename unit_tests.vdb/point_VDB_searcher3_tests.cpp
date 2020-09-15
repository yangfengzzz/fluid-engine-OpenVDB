//
//  point_VDB_searcher3_tests.cpp
//  unittest
//
//  Created by Feng Yang on 2020/1/5.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.vdb/vdb_point_searcher3.hpp"

#include "../external/gtest/include/gtest/gtest.h"
#include <vector>

using namespace vdb;

TEST(PointVDBSearcher3, ForEachNearbyPoint) {
    vox::Array1<vox::Vector3D> points = {
        vox::Vector3D(0, 1, 3),
        vox::Vector3D(2, 5, 4),
        vox::Vector3D(-1, 3, 0)
    };
    
    PointVDBSearcher3 searcher;
    searcher.build(points);
    
    int cnt = 0;
    searcher.forEachNearbyPoint(
                                vox::Vector3D(0, 0, 0),
                                std::sqrt(10.0),
                                [&](size_t i, const vox::Vector3D& pt) {
        EXPECT_TRUE(i == 0 || i == 2);
        
        if (i == 0) {
            EXPECT_EQ(points[0], pt);
        } else if (i == 2) {
            EXPECT_EQ(points[2], pt);
        }
        
        ++cnt;
    });
    
    EXPECT_EQ(2, cnt);
}

TEST(PointVDBSearcher3, ForEachNearbyPointEmpty) {
    vox::Array1<vox::Vector3D> points;
    
    PointVDBSearcher3 searcher;
    searcher.build(points.accessor());
    
    searcher.forEachNearbyPoint(vox::Vector3D(0, 0, 0), std::sqrt(10.0),
                                [](size_t, const vox::Vector3D&) {});
}

TEST(PointVDBSearcher3, CopyConstructor) {
    vox::Array1<vox::Vector3D> points = {
        vox::Vector3D(0, 1, 3),
        vox::Vector3D(2, 5, 4),
        vox::Vector3D(-1, 3, 0)
    };
    
    PointVDBSearcher3 searcher;
    searcher.build(points);
    
    PointVDBSearcher3 searcher2(searcher);
    int cnt = 0;
    searcher2.forEachNearbyPoint(
                                 vox::Vector3D(0, 0, 0),
                                 std::sqrt(10.0),
                                 [&](size_t i, const vox::Vector3D& pt) {
        EXPECT_TRUE(i == 0 || i == 2);
        
        if (i == 0) {
            EXPECT_EQ(points[0], pt);
        } else if (i == 2) {
            EXPECT_EQ(points[2], pt);
        }
        
        ++cnt;
    });
    
    EXPECT_EQ(2, cnt);
}

