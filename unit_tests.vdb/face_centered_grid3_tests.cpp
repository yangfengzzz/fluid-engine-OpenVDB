// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.vdb/vdb_face_centered_grid3.h"
#include "../src.vdb/vdb_math_utils.h"
#include "../src.vdb/vdb_helper.h"
#include "../external/gtest/include/gtest/gtest.h"
#include <vector>

using namespace vdb;

TEST(FaceCenteredGrid3, Constructors) {
    // Default constructors
    FaceCenteredGrid3 grid1;
    EXPECT_EQ(0u, grid1.resolution().x);
    EXPECT_EQ(0u, grid1.resolution().y);
    EXPECT_EQ(0u, grid1.resolution().z);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().x);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().y);
    EXPECT_DOUBLE_EQ(1.0, grid1.gridSpacing().z);
    EXPECT_DOUBLE_EQ(0.0, grid1.origin().x);
    EXPECT_DOUBLE_EQ(0.0, grid1.origin().y);
    EXPECT_DOUBLE_EQ(0.0, grid1.origin().z);
    EXPECT_EQ(0u, grid1.uSize().x);
    EXPECT_EQ(0u, grid1.uSize().y);
    EXPECT_EQ(0u, grid1.uSize().z);
    EXPECT_EQ(0u, grid1.vSize().x);
    EXPECT_EQ(0u, grid1.vSize().y);
    EXPECT_EQ(0u, grid1.vSize().z);
    EXPECT_EQ(0u, grid1.wSize().x);
    EXPECT_EQ(0u, grid1.wSize().y);
    EXPECT_EQ(0u, grid1.wSize().z);
    EXPECT_DOUBLE_EQ(0.0, grid1.uOrigin().x);
    EXPECT_DOUBLE_EQ(0.5, grid1.uOrigin().y);
    EXPECT_DOUBLE_EQ(0.5, grid1.uOrigin().z);
    EXPECT_DOUBLE_EQ(0.5, grid1.vOrigin().x);
    EXPECT_DOUBLE_EQ(0.0, grid1.vOrigin().y);
    EXPECT_DOUBLE_EQ(0.5, grid1.vOrigin().z);
    EXPECT_DOUBLE_EQ(0.5, grid1.wOrigin().x);
    EXPECT_DOUBLE_EQ(0.5, grid1.wOrigin().y);
    EXPECT_DOUBLE_EQ(0.0, grid1.wOrigin().z);
    
    // Constructor with params
    FaceCenteredGrid3 grid2(
                            5, 4, 3,
                            1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0,
                            7.0, 8.0, 9.0);
    EXPECT_EQ(5u, grid2.resolution().x);
    EXPECT_EQ(4u, grid2.resolution().y);
    EXPECT_EQ(3u, grid2.resolution().z);
    EXPECT_DOUBLE_EQ(1.0, grid2.gridSpacing().x);
    EXPECT_DOUBLE_EQ(2.0, grid2.gridSpacing().y);
    EXPECT_DOUBLE_EQ(3.0, grid2.gridSpacing().z);
    EXPECT_DOUBLE_EQ(4.0, grid2.origin().x);
    EXPECT_DOUBLE_EQ(5.0, grid2.origin().y);
    EXPECT_DOUBLE_EQ(6.0, grid2.origin().z);
    EXPECT_EQ(6u, grid2.uSize().x);
    EXPECT_EQ(4u, grid2.uSize().y);
    EXPECT_EQ(3u, grid2.uSize().z);
    EXPECT_EQ(5u, grid2.vSize().x);
    EXPECT_EQ(5u, grid2.vSize().y);
    EXPECT_EQ(3u, grid2.vSize().z);
    EXPECT_EQ(5u, grid2.wSize().x);
    EXPECT_EQ(4u, grid2.wSize().y);
    EXPECT_EQ(4u, grid2.wSize().z);
    EXPECT_DOUBLE_EQ(4.0, grid2.uOrigin().x);
    EXPECT_DOUBLE_EQ(6.0, grid2.uOrigin().y);
    EXPECT_DOUBLE_EQ(7.5, grid2.uOrigin().z);
    EXPECT_DOUBLE_EQ(4.5, grid2.vOrigin().x);
    EXPECT_DOUBLE_EQ(5.0, grid2.vOrigin().y);
    EXPECT_DOUBLE_EQ(7.5, grid2.vOrigin().z);
    EXPECT_DOUBLE_EQ(4.5, grid2.wOrigin().x);
    EXPECT_DOUBLE_EQ(6.0, grid2.wOrigin().y);
    EXPECT_DOUBLE_EQ(6.0, grid2.wOrigin().z);
    grid2.forEachUIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(7.0, grid2.u(coord));
    });
    grid2.forEachVIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(8.0, grid2.v(coord));
    });
    grid2.forEachWIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(9.0, grid2.w(coord));
    });
    
    // Copy constructor
    FaceCenteredGrid3 grid3(grid2);
    EXPECT_EQ(5u, grid3.resolution().x);
    EXPECT_EQ(4u, grid3.resolution().y);
    EXPECT_EQ(3u, grid3.resolution().z);
    EXPECT_DOUBLE_EQ(1.0, grid3.gridSpacing().x);
    EXPECT_DOUBLE_EQ(2.0, grid3.gridSpacing().y);
    EXPECT_DOUBLE_EQ(3.0, grid3.gridSpacing().z);
    EXPECT_DOUBLE_EQ(4.0, grid3.origin().x);
    EXPECT_DOUBLE_EQ(5.0, grid3.origin().y);
    EXPECT_DOUBLE_EQ(6.0, grid3.origin().z);
    EXPECT_EQ(6u, grid3.uSize().x);
    EXPECT_EQ(4u, grid3.uSize().y);
    EXPECT_EQ(3u, grid3.uSize().z);
    EXPECT_EQ(5u, grid3.vSize().x);
    EXPECT_EQ(5u, grid3.vSize().y);
    EXPECT_EQ(3u, grid3.vSize().z);
    EXPECT_EQ(5u, grid3.wSize().x);
    EXPECT_EQ(4u, grid3.wSize().y);
    EXPECT_EQ(4u, grid3.wSize().z);
    EXPECT_DOUBLE_EQ(4.0, grid3.uOrigin().x);
    EXPECT_DOUBLE_EQ(6.0, grid3.uOrigin().y);
    EXPECT_DOUBLE_EQ(7.5, grid3.uOrigin().z);
    EXPECT_DOUBLE_EQ(4.5, grid3.vOrigin().x);
    EXPECT_DOUBLE_EQ(5.0, grid3.vOrigin().y);
    EXPECT_DOUBLE_EQ(7.5, grid3.vOrigin().z);
    EXPECT_DOUBLE_EQ(4.5, grid3.wOrigin().x);
    EXPECT_DOUBLE_EQ(6.0, grid3.wOrigin().y);
    EXPECT_DOUBLE_EQ(6.0, grid3.wOrigin().z);
    grid3.forEachUIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(7.0, grid3.u(coord));
    });
    grid3.forEachVIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(8.0, grid3.v(coord));
    });
    grid3.forEachWIndex([&] (const openvdb::Coord& coord) {
        EXPECT_DOUBLE_EQ(9.0, grid3.w(coord));
    });
}

TEST(FaceCenteredGrid3, Fill) {
    FaceCenteredGrid3 grid(
                           5, 4, 6,
                           1.0, 1.0, 1.0,
                           0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0);
    grid.fill(vox::Vector3D(42.0, 27.0, 31.0),
              vox::ExecutionPolicy::kSerial);
    
    for (uint k = 0; k < grid.uSize().z; ++k) {
        for (uint j = 0; j < grid.uSize().y; ++j) {
            for (uint i = 0; i < grid.uSize().x; ++i) {
                EXPECT_DOUBLE_EQ(42.0, grid.u(openvdb::Coord(i, j, k)));
            }
        }
    }
    
    for (uint k = 0; k < grid.vSize().z; ++k) {
        for (uint j = 0; j < grid.vSize().y; ++j) {
            for (uint i = 0; i < grid.vSize().x; ++i) {
                EXPECT_DOUBLE_EQ(27.0, grid.v(openvdb::Coord(i, j, k)));
            }
        }
    }
    
    for (uint k = 0; k < grid.wSize().z; ++k) {
        for (uint j = 0; j < grid.wSize().y; ++j) {
            for (uint i = 0; i < grid.wSize().x; ++i) {
                EXPECT_DOUBLE_EQ(31.0, grid.w(openvdb::Coord(i, j, k)));
            }
        }
    }
    
    auto func = [](const vox::Vector3D& x) { return x; };
    grid.fill(func, vox::ExecutionPolicy::kSerial);
    
    for (uint k = 0; k < grid.uSize().z; ++k) {
        for (uint j = 0; j < grid.uSize().y; ++j) {
            for (uint i = 0; i < grid.uSize().x; ++i) {
                EXPECT_DOUBLE_EQ(static_cast<double>(i),
                                 grid.u(openvdb::Coord(i, j, k)));
            }
        }
    }
    
    for (uint k = 0; k < grid.vSize().z; ++k) {
        for (uint j = 0; j < grid.vSize().y; ++j) {
            for (uint i = 0; i < grid.vSize().x; ++i) {
                EXPECT_DOUBLE_EQ(static_cast<double>(j),
                                 grid.v(openvdb::Coord(i, j, k)));
            }
        }
    }
    
    for (uint k = 0; k < grid.wSize().z; ++k) {
        for (uint j = 0; j < grid.wSize().y; ++j) {
            for (uint i = 0; i < grid.wSize().x; ++i) {
                EXPECT_DOUBLE_EQ(static_cast<double>(k),
                                 grid.w(openvdb::Coord(i, j, k)));
            }
        }
    }
}

TEST(FaceCenteredGrid3, DivergenceAtCellCenter) {
    FaceCenteredGrid3 grid(5, 8, 6);
    
    grid.fill(vox::Vector3D(1.0, -2.0, 3.0),
              vox::ExecutionPolicy::kSerial);
    
    for (uint k = 0; k < grid.resolution().z; ++k) {
        for (uint j = 0; j < grid.resolution().y; ++j) {
            for (uint i = 0; i < grid.resolution().x; ++i) {
                EXPECT_DOUBLE_EQ(0.0, grid.divergenceAtCellCenter(i, j, k));
            }
        }
    }
    
    grid.fill([](const vox::Vector3D& x) { return x; },
              vox::ExecutionPolicy::kSerial);
    
    for (uint k = 0; k < grid.resolution().z; ++k) {
        for (uint j = 0; j < grid.resolution().y; ++j) {
            for (uint i = 0; i < grid.resolution().x; ++i) {
                EXPECT_DOUBLE_EQ(3.0, grid.divergenceAtCellCenter(i, j, k));
            }
        }
    }
}

TEST(FaceCenteredGrid3, CurlAtCellCenter) {
    FaceCenteredGrid3 grid(5, 8, 6, 2.0, 3.0, 1.5);
    
    grid.fill(vox::Vector3D(1.0, -2.0, 3.0),
              vox::ExecutionPolicy::kSerial);
    
    for (uint k = 0; k < grid.resolution().z; ++k) {
        for (uint j = 0; j < grid.resolution().y; ++j) {
            for (uint i = 0; i < grid.resolution().x; ++i) {
                vox::Vector3D curl = grid.curlAtCellCenter(i, j, k);
                EXPECT_DOUBLE_EQ(0.0, curl.x);
                EXPECT_DOUBLE_EQ(0.0, curl.y);
                EXPECT_DOUBLE_EQ(0.0, curl.z);
            }
        }
    }
    
    grid.fill([](const vox::Vector3D& x) { return vox::Vector3D(x.y, x.z, x.x); },
              vox::ExecutionPolicy::kSerial);
    
    for (uint k = 1; k < grid.resolution().z - 1; ++k) {
        for (uint j = 1; j < grid.resolution().y - 1; ++j) {
            for (uint i = 1; i < grid.resolution().x - 1; ++i) {
                vox::Vector3D curl = grid.curlAtCellCenter(i, j, k);
                EXPECT_DOUBLE_EQ(-1.0, curl.x);
                EXPECT_DOUBLE_EQ(-1.0, curl.y);
                EXPECT_DOUBLE_EQ(-1.0, curl.z);
            }
        }
    }
}

TEST(FaceCenteredGrid3, ValueAtCellCenter) {
    FaceCenteredGrid3 grid(5, 8, 6, 2.0, 3.0, 1.5);
    grid.fill([&](const vox::Vector3D& x) {
        return vox::Vector3D(3.0 * x.y + 1.0, 5.0 * x.z + 7.0, -1.0 * x.x - 9.0);
    }, vox::ExecutionPolicy::kSerial);
    
    auto pos = grid.cellCenterPosition();
    grid.forEachCellIndex([&](uint i, uint j, uint k) {
        vox::Vector3D val = grid.valueAtCellCenter(i, j, k);
        vox::Vector3D x = pos(openvdb::Coord(i, j, k));
        vox::Vector3D expected
        = vox::Vector3D(3.0 * x.y + 1.0, 5.0 * x.z + 7.0, -1.0 * x.x - 9.0);
        EXPECT_NEAR(expected.x, val.x, 1e-6);
        EXPECT_NEAR(expected.y, val.y, 1e-6);
        EXPECT_NEAR(expected.z, val.z, 1e-6);
    });
}

TEST(FaceCenteredGrid3, Sample) {
    FaceCenteredGrid3 grid(5, 8, 6, 2.0, 3.0, 1.5);
    grid.fill([&](const vox::Vector3D& x) {
        return vox::Vector3D(3.0 * x.y + 1.0, 5.0 * x.z + 7.0, -1.0 * x.x - 9.0);
    }, vox::ExecutionPolicy::kSerial);
    
    auto pos = grid.cellCenterPosition();
    grid.forEachCellIndex([&](uint i, uint j, uint k) {
        vox::Vector3D x = pos(openvdb::Coord(i, j, k));
        vox::Vector3D val = grid.sample(x);
        vox::Vector3D expected
        = vox::Vector3D(3.0 * x.y + 1.0, 5.0 * x.z + 7.0, -1.0 * x.x - 9.0);
        EXPECT_NEAR(expected.x, val.x, 1e-6);
        EXPECT_NEAR(expected.y, val.y, 1e-6);
        EXPECT_NEAR(expected.z, val.z, 1e-6);
    });
}

TEST(FaceCenteredGrid3, Builder) {
    {
        auto builder = FaceCenteredGrid3::builder();
        
        auto grid = builder.build(
                                  vox::Size3(5, 2, 7),
                                  vox::Vector3D(2.0, 4.0, 1.5),
                                  vox::Vector3D(-1.0, 2.0, 7.0),
                                  vox::Vector3D(3.0, 5.0, -2.0));
        EXPECT_EQ(vox::Size3(5, 2, 7), grid->resolution());
        EXPECT_EQ(vox::Vector3D(2.0, 4.0, 1.5), grid->gridSpacing());
        EXPECT_EQ(vox::Vector3D(-1.0, 2.0, 7.0), grid->origin());
        
        auto faceCenteredGrid
        = std::dynamic_pointer_cast<FaceCenteredGrid3>(grid);
        faceCenteredGrid->forEachUIndex(
                                        [&faceCenteredGrid](const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(3.0, faceCenteredGrid->u(coord));
        });
        faceCenteredGrid->forEachVIndex(
                                        [&faceCenteredGrid](const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(5.0, faceCenteredGrid->v(coord));
        });
        faceCenteredGrid->forEachWIndex(
                                        [&faceCenteredGrid](const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(-2.0, faceCenteredGrid->w(coord));
        });
    }
    
    {
        auto grid = FaceCenteredGrid3::builder()
        .withResolution(5, 2, 7)
        .withGridSpacing(2.0, 4.0, 1.5)
        .withOrigin(-1.0, 2.0, 7.0)
        .withInitialValue(3.0, 5.0, -2.0)
        .build();
        EXPECT_EQ(vox::Size3(5, 2, 7), grid.resolution());
        EXPECT_EQ(vox::Vector3D(2.0, 4.0, 1.5), grid.gridSpacing());
        EXPECT_EQ(vox::Vector3D(-1.0, 2.0, 7.0), grid.origin());
        
        grid.forEachUIndex(
                           [&](const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(3.0, grid.u(coord));
        });
        grid.forEachVIndex(
                           [&](const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(5.0, grid.v(coord));
        });
        grid.forEachWIndex(
                           [&](const openvdb::Coord& coord) {
            EXPECT_DOUBLE_EQ(-2.0, grid.w(coord));
        });
    }
}
