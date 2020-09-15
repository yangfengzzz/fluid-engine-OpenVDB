//
//  vdb_emitter3_tests.cpp
//  vdb_tests
//
//  Created by Feng Yang on 2020/2/6.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_cell_centered_vector_grid3.h"
#include "../src.common/custom_vector_field3.h"
#include "../src.common/level_set_utils.h"
#include "../src.common/sphere3.h"
#include "../src.common/triangle_mesh3.h"
#include "../src.common/implicit_triangle_mesh3.h"
#include "../src.common/surface_to_implicit3.h"
#include "../src.vdb/vdb_volume_grid_emitter3.hpp"

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(VolumeGridEmitter3, Velocity) {
    auto sphere = vox::Sphere3::builder()
    .withCenter({0.5, 0.75, 0.5})
    .withRadius(0.15)
    .makeShared();
    
    auto emitter = GridEmitter3::builder()
    .withSourceRegion(sphere)
    .makeShared();
    
    auto grid = CellCenteredVectorGrid3::builder()
    .withResolution({16, 16, 16})
    .withGridSpacing({1.0/16.0, 1.0/16.0, 1.0/16.0})
    .withOrigin({0, 0, 0})
    .makeShared();
    grid->getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                  openvdb::Coord(16, 16, 16)),
                               openvdb::Vec3d::zero(), true);
    
    auto mapper = [] (double sdf, const vox::Vector3D& pt,
                      const vox::Vector3D& oldVal) {
        if (sdf < 0.0) {
            return vox::Vector3D(pt.y, -pt.x, 3.5);
        } else {
            return vox::Vector3D(oldVal);
        }
    };
    
    emitter->addTarget(grid, mapper);
    
    emitter->update(0.0, 0.01);
    
    auto pos = grid->dataPosition();
    grid->forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        vox::Vector3D gx = pos(coord);
        double sdf = emitter->sourceRegion()->signedDistance(gx);
        if (vox::isInsideSdf(sdf)) {
            openvdb::Vec3d answer{gx.y, -gx.x, 3.5};
            openvdb::Vec3d acttual = (*grid)(coord);
            
            EXPECT_NEAR(answer.x(), acttual.x(), 1e-6);
            EXPECT_NEAR(answer.y(), acttual.y(), 1e-6);
        }
    });
}

TEST(VolumeGridEmitter3, SignedDistance) {
    auto sphere = vox::Sphere3::builder()
    .withCenter({0.5, 0.75, 0.5})
    .withRadius(0.15)
    .makeShared();
    
    auto emitter = GridEmitter3::builder()
    .withSourceRegion(sphere)
    .makeShared();
    
    auto grid = CellCenteredScalarGrid3::builder()
    .withResolution({16, 16, 16})
    .withGridSpacing({1.0/16.0, 1.0/16.0, 1.0/16.0})
    .withOrigin({0, 0, 0})
    .withInitialValue(vox::kMaxD)
    .makeShared();
    grid->getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                  openvdb::Coord(16, 16, 16)),
                               vox::kMaxD, true);
    
    emitter->addSignedDistanceTarget(grid);
    
    emitter->update(0.0, 0.01);
    
    auto pos = grid->dataPosition();
    grid->forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        vox::Vector3D gx = pos(coord);
        double answer = (sphere->center - gx).length() - 0.15;
        double acttual = (*grid)(coord);
        
        EXPECT_NEAR(answer, acttual, 1e-6);
    });
}

TEST(VolumeGridEmitter3, StepFunction) {
    auto sphere = vox::Sphere3::builder()
    .withCenter({0.5, 0.75, 0.5})
    .withRadius(0.15)
    .makeShared();
    
    auto emitter = GridEmitter3::builder()
    .withSourceRegion(sphere)
    .makeShared();
    
    auto grid = CellCenteredScalarGrid3::builder()
    .withResolution({16, 16, 16})
    .withGridSpacing({1.0/16.0, 1.0/16.0, 1.0/16.0})
    .withOrigin({0, 0, 0})
    .makeShared();
    grid->getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                  openvdb::Coord(16, 16, 16)),
                               0.0, true);
    
    emitter->addStepFunctionTarget(grid, 3.0, 7.0);
    
    emitter->update(0.0, 0.01);
    
    auto pos = grid->dataPosition();
    grid->forEachDataPointIndex([&] (const openvdb::Coord& coord) {
        vox::Vector3D gx = pos(coord);
        double answer = (sphere->center - gx).length() - 0.15;
        answer = 4.0 * (1.0 - vox::smearedHeavisideSdf(answer * 16.0)) + 3.0;
        double acttual = (*grid)(coord);
        
        EXPECT_NEAR(answer, acttual, 1e-6);
    });
}
