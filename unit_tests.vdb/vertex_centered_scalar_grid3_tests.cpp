// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.vdb/vdb_vertex_centered_scalar_grid3.h"
#include "../src.vdb/vdb_math_utils.h"
#include "../src.vdb/vdb_helper.h"
#include "../external/gtest/include/gtest/gtest.h"
#include <vector>

using namespace vdb;

TEST(VertexCenteredScalarGrid3, Constructors) {
    // Default constructors
    VertexCenteredScalarGrid3 grid1;
    EXPECT_EQ(0u, grid1.resolution().x);
    EXPECT_EQ(0u, grid1.resolution().y);
    EXPECT_EQ(0u, grid1.resolution().z);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().x);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().y);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().z);
    EXPECT_DOUBLE_EQ(0.0, grid1.origin().x);
    EXPECT_DOUBLE_EQ(0.0, grid1.origin().y);
    EXPECT_DOUBLE_EQ(0.0, grid1.origin().z);
    EXPECT_EQ(0u, grid1.dataSize().x);
    EXPECT_EQ(0u, grid1.dataSize().y);
    EXPECT_EQ(0u, grid1.dataSize().z);
    EXPECT_DOUBLE_EQ(0.0, grid1.dataOrigin().x);
    EXPECT_DOUBLE_EQ(0.0, grid1.dataOrigin().y);
    EXPECT_DOUBLE_EQ(0.0, grid1.dataOrigin().z);
    
    // Constructor with params
    VertexCenteredScalarGrid3 grid2(5, 4, 3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
    EXPECT_EQ(5u, grid2.resolution().x);
    EXPECT_EQ(4u, grid2.resolution().y);
    EXPECT_EQ(3u, grid2.resolution().z);
    EXPECT_DOUBLE_EQ(1.0, grid2.gridSpacing().x);
    EXPECT_DOUBLE_EQ(2.0, grid2.gridSpacing().y);
    EXPECT_DOUBLE_EQ(3.0, grid2.gridSpacing().z);
    EXPECT_DOUBLE_EQ(4.0, grid2.origin().x);
    EXPECT_DOUBLE_EQ(5.0, grid2.origin().y);
    EXPECT_DOUBLE_EQ(6.0, grid2.origin().z);
    EXPECT_EQ(6u, grid2.dataSize().x);
    EXPECT_EQ(5u, grid2.dataSize().y);
    EXPECT_EQ(4u, grid2.dataSize().z);
    EXPECT_DOUBLE_EQ(4.0, grid2.dataOrigin().x);
    EXPECT_DOUBLE_EQ(5.0, grid2.dataOrigin().y);
    EXPECT_DOUBLE_EQ(6.0, grid2.dataOrigin().z);
    grid2.forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(7.0, grid2(coord));
    });
    
    // Copy constructor
    VertexCenteredScalarGrid3 grid3(grid2);
    EXPECT_EQ(5u, grid3.resolution().x);
    EXPECT_EQ(4u, grid3.resolution().y);
    EXPECT_EQ(3u, grid3.resolution().z);
    EXPECT_DOUBLE_EQ(1.0, grid3.gridSpacing().x);
    EXPECT_DOUBLE_EQ(2.0, grid3.gridSpacing().y);
    EXPECT_DOUBLE_EQ(3.0, grid3.gridSpacing().z);
    EXPECT_DOUBLE_EQ(4.0, grid3.origin().x);
    EXPECT_DOUBLE_EQ(5.0, grid3.origin().y);
    EXPECT_DOUBLE_EQ(6.0, grid3.origin().z);
    EXPECT_EQ(6u, grid3.dataSize().x);
    EXPECT_EQ(5u, grid3.dataSize().y);
    EXPECT_EQ(4u, grid3.dataSize().z);
    EXPECT_DOUBLE_EQ(4.0, grid3.dataOrigin().x);
    EXPECT_DOUBLE_EQ(5.0, grid3.dataOrigin().y);
    EXPECT_DOUBLE_EQ(6.0, grid3.dataOrigin().z);
    grid3.forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(7.0, grid3(coord));
    });
}

TEST(VertexCenteredScalarGrid3, Swap) {
    VertexCenteredScalarGrid3 grid1(5, 4, 3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
    VertexCenteredScalarGrid3 grid2(3, 8, 5, 2.0, 3.0, 1.0, 5.0, 4.0, 7.0, 8.0);
    grid1.swap(&grid2);
    
    EXPECT_EQ(3u, grid1.resolution().x);
    EXPECT_EQ(8u, grid1.resolution().y);
    EXPECT_EQ(5u, grid1.resolution().z);
    EXPECT_DOUBLE_EQ(2.0, grid1.gridSpacing().x);
    EXPECT_DOUBLE_EQ(3.0, grid1.gridSpacing().y);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().z);
    EXPECT_DOUBLE_EQ(5.0, grid1.origin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1.origin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1.origin().z);
    EXPECT_EQ(4u, grid1.dataSize().x);
    EXPECT_EQ(9u, grid1.dataSize().y);
    EXPECT_EQ(6u, grid1.dataSize().z);
    EXPECT_DOUBLE_EQ(5.0, grid1.dataOrigin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1.dataOrigin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1.dataOrigin().z);
    grid1.forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(8.0, grid1(coord));
    });
    
    EXPECT_EQ(5u, grid2.resolution().x);
    EXPECT_EQ(4u, grid2.resolution().y);
    EXPECT_EQ(3u, grid2.resolution().z);
    EXPECT_DOUBLE_EQ(1.0, grid2.gridSpacing().x);
    EXPECT_DOUBLE_EQ(2.0, grid2.gridSpacing().y);
    EXPECT_DOUBLE_EQ(3.0, grid2.gridSpacing().z);
    EXPECT_DOUBLE_EQ(4.0, grid2.origin().x);
    EXPECT_DOUBLE_EQ(5.0, grid2.origin().y);
    EXPECT_DOUBLE_EQ(6.0, grid2.origin().z);
    EXPECT_EQ(6u, grid2.dataSize().x);
    EXPECT_EQ(5u, grid2.dataSize().y);
    EXPECT_EQ(4u, grid2.dataSize().z);
    EXPECT_DOUBLE_EQ(4.0, grid2.dataOrigin().x);
    EXPECT_DOUBLE_EQ(5.0, grid2.dataOrigin().y);
    EXPECT_DOUBLE_EQ(6.0, grid2.dataOrigin().z);
    grid2.forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(7.0, grid2(coord));
    });
}

TEST(VertexCenteredScalarGrid3, Set) {
    VertexCenteredScalarGrid3 grid1(5, 4, 3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
    VertexCenteredScalarGrid3 grid2(3, 8, 5, 2.0, 3.0, 1.0, 5.0, 4.0, 7.0, 8.0);
    grid1.set(grid2);
    
    EXPECT_EQ(3u, grid1.resolution().x);
    EXPECT_EQ(8u, grid1.resolution().y);
    EXPECT_EQ(5u, grid1.resolution().z);
    EXPECT_DOUBLE_EQ(2.0, grid1.gridSpacing().x);
    EXPECT_DOUBLE_EQ(3.0, grid1.gridSpacing().y);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().z);
    EXPECT_DOUBLE_EQ(5.0, grid1.origin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1.origin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1.origin().z);
    EXPECT_EQ(4u, grid1.dataSize().x);
    EXPECT_EQ(9u, grid1.dataSize().y);
    EXPECT_EQ(6u, grid1.dataSize().z);
    EXPECT_DOUBLE_EQ(5.0, grid1.dataOrigin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1.dataOrigin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1.dataOrigin().z);
    grid1.forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(8.0, grid1(coord));
    });
}

TEST(VertexCenteredScalarGrid3, AssignmentOperator) {
    VertexCenteredScalarGrid3 grid1(5, 4, 3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
    VertexCenteredScalarGrid3 grid2(3, 8, 5, 2.0, 3.0, 1.0, 5.0, 4.0, 7.0, 8.0);
    grid1 = grid2;
    
    EXPECT_EQ(3u, grid1.resolution().x);
    EXPECT_EQ(8u, grid1.resolution().y);
    EXPECT_EQ(5u, grid1.resolution().z);
    EXPECT_DOUBLE_EQ(2.0, grid1.gridSpacing().x);
    EXPECT_DOUBLE_EQ(3.0, grid1.gridSpacing().y);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().z);
    EXPECT_DOUBLE_EQ(5.0, grid1.origin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1.origin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1.origin().z);
    EXPECT_EQ(4u, grid1.dataSize().x);
    EXPECT_EQ(9u, grid1.dataSize().y);
    EXPECT_EQ(6u, grid1.dataSize().z);
    EXPECT_DOUBLE_EQ(5.0, grid1.dataOrigin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1.dataOrigin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1.dataOrigin().z);
    grid1.forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(8.0, grid1(coord));
    });
}

TEST(VertexCenteredScalarGrid3, Clone) {
    VertexCenteredScalarGrid3 grid2(3, 8, 5, 2.0, 3.0, 1.0, 5.0, 4.0, 7.0, 8.0);
    auto grid1 = grid2.clone();
    
    EXPECT_EQ(3u, grid1->resolution().x);
    EXPECT_EQ(8u, grid1->resolution().y);
    EXPECT_EQ(5u, grid1->resolution().z);
    EXPECT_DOUBLE_EQ(2.0, grid1->gridSpacing().x);
    EXPECT_DOUBLE_EQ(3.0, grid1->gridSpacing().y);
    EXPECT_DOUBLE_EQ(1.0, grid1->gridSpacing().z);
    EXPECT_DOUBLE_EQ(5.0, grid1->origin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1->origin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1->origin().z);
    EXPECT_EQ(4u, grid1->dataSize().x);
    EXPECT_EQ(9u, grid1->dataSize().y);
    EXPECT_EQ(6u, grid1->dataSize().z);
    EXPECT_DOUBLE_EQ(5.0, grid1->dataOrigin().x);
    EXPECT_DOUBLE_EQ(4.0, grid1->dataOrigin().y);
    EXPECT_DOUBLE_EQ(7.0, grid1->dataOrigin().z);
    grid1->forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(8.0, (*grid1)(coord));
    });
}

TEST(VertexCenteredScalarGrid3, Builder) {
    {
        auto grid1 = VertexCenteredScalarGrid3::builder().build(
                                                                {3, 8, 5}, {2.0, 3.0, 1.0}, {5.0, 4.0, 7.0}, 8.0);
        
        auto grid2
        = std::dynamic_pointer_cast<VertexCenteredScalarGrid3>(grid1);
        EXPECT_TRUE(grid2 != nullptr);
        
        EXPECT_EQ(3u, grid1->resolution().x);
        EXPECT_EQ(8u, grid1->resolution().y);
        EXPECT_EQ(5u, grid1->resolution().z);
        EXPECT_DOUBLE_EQ(2.0, grid1->gridSpacing().x);
        EXPECT_DOUBLE_EQ(3.0, grid1->gridSpacing().y);
        EXPECT_DOUBLE_EQ(1.0, grid1->gridSpacing().z);
        EXPECT_DOUBLE_EQ(5.0, grid1->origin().x);
        EXPECT_DOUBLE_EQ(4.0, grid1->origin().y);
        EXPECT_DOUBLE_EQ(7.0, grid1->origin().z);
        EXPECT_EQ(4u, grid1->dataSize().x);
        EXPECT_EQ(9u, grid1->dataSize().y);
        EXPECT_EQ(6u, grid1->dataSize().z);
        EXPECT_DOUBLE_EQ(5.0, grid1->dataOrigin().x);
        EXPECT_DOUBLE_EQ(4.0, grid1->dataOrigin().y);
        EXPECT_DOUBLE_EQ(7.0, grid1->dataOrigin().z);
        grid1->forEachDataPointIndex([&] (const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(8.0, (*grid1)(coord));
        });
    }
    
    {
        auto grid1 = VertexCenteredScalarGrid3::builder()
        .withResolution(3, 8, 5)
        .withGridSpacing(2, 3, 1)
        .withOrigin(5, 4, 7)
        .withInitialValue(8.0)
        .build();
        
        EXPECT_EQ(3u, grid1.resolution().x);
        EXPECT_EQ(8u, grid1.resolution().y);
        EXPECT_EQ(5u, grid1.resolution().z);
        EXPECT_DOUBLE_EQ(2.0, grid1.gridSpacing().x);
        EXPECT_DOUBLE_EQ(3.0, grid1.gridSpacing().y);
        EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().z);
        EXPECT_DOUBLE_EQ(5.0, grid1.origin().x);
        EXPECT_DOUBLE_EQ(4.0, grid1.origin().y);
        EXPECT_DOUBLE_EQ(7.0, grid1.origin().z);
        EXPECT_EQ(4u, grid1.dataSize().x);
        EXPECT_EQ(9u, grid1.dataSize().y);
        EXPECT_EQ(6u, grid1.dataSize().z);
        EXPECT_DOUBLE_EQ(5.0, grid1.dataOrigin().x);
        EXPECT_DOUBLE_EQ(4.0, grid1.dataOrigin().y);
        EXPECT_DOUBLE_EQ(7.0, grid1.dataOrigin().z);
        grid1.forEachDataPointIndex([&] (const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(8.0, grid1(coord));
        });
    }
}

TEST(VertexCenteredScalarGrid3, Fill) {
    VertexCenteredScalarGrid3 grid(5, 4, 6, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    grid.fill(42.0, vox::ExecutionPolicy::kSerial);
    
    for (uint k = 0; k < grid.dataSize().z; ++k) {
        for (uint j = 0; j < grid.dataSize().y; ++j) {
            for (uint i = 0; i < grid.dataSize().x; ++i) {
                EXPECT_DOUBLE_EQ(42.0, grid(openvdb::Coord(i, j, k)));
            }
        }
    }
    
    auto func = [](const vox::Vector3D& x) {
        if (x.x < 3.0) {
            return 2.0;
        } else {
            return 5.0;
        }
    };
    
    grid.fill(func, vox::ExecutionPolicy::kSerial);
    
    for (uint k = 0; k < grid.dataSize().z; ++k) {
        for (uint j = 0; j < grid.dataSize().y; ++j) {
            for (uint i = 0; i < grid.dataSize().x; ++i) {
                if (i < 3) {
                    EXPECT_DOUBLE_EQ(2.0, grid(openvdb::Coord(i, j, k)));
                } else {
                    EXPECT_DOUBLE_EQ(5.0, grid(openvdb::Coord(i, j, k)));
                }
            }
        }
    }
}
