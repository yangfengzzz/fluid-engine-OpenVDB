//
//  advection_solvers_tests.cpp
//  vdb_manual_tests
//
//  Created by Feng Yang on 2020/2/12.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "manual_tests.h"
#include "../src.vdb/vdb_helper.h"
#include "../src.common/array2.h"
#include "../src.common/array3.h"
#include "../src.common/box3.h"
#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_cell_centered_vector_grid3.h"
#include "../src.common/constant_vector_field3.h"
#include "../src.common/constants.h"
#include "../src.vdb/vdb_cubic_semi_lagrangian3.hpp"
#include "../src.common/custom_scalar_field3.h"
#include "../src.common/custom_vector_field3.h"
#include "../src.vdb/vdb_semi_lagrangian3.hpp"
#include "../src.common/triangle_mesh3.h"

#include <algorithm>

using namespace vdb;

JET_TESTS(SemiLagrangian3);

JET_BEGIN_TEST_F(SemiLagrangian3, Boundary) {
    CellCenteredVectorGrid3 src(200, 200, 200,
                                1.0/200.0, 1.0/200.0, 1.0/200.0);
    CellCenteredVectorGrid3 dst(200, 200, 200,
                                1.0/200.0, 1.0/200.0, 1.0/200.0);
    src.getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                openvdb::Coord(199, 199, 199)),
                             openvdb::Vec3d::zero(), true);
    dst.getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                openvdb::Coord(199, 199, 199)),
                             openvdb::Vec3d::zero(), true);
    
    src.fill([&](const vox::Vector3D& pt) -> vox::Vector3D {
        return {
            0.5 * (std::sin(15 * pt.x) + 1.0),
            0.5 * (std::sin(15 * pt.y) + 1.0),
            0.0};
    }, vox::ExecutionPolicy::kSerial);

    vox::ConstantVectorField3 flow(vox::Vector3D(1.0, 1.0, 0.0));
    vox::CustomScalarField3 boundarySdf([](const vox::Vector3D& pt) {
        return vox::Vector3D(0.5, 0.5, 0.5).distanceTo(pt) - 0.25;
    });
    
    vox::Array3<double> data(3, src.resolution().x, src.resolution().y);
    data.forEachIndex([&](uint i, uint j, uint k) {
        if (i < 2) {
            data(i, j, k) = src.getGrid()->tree().getValue(openvdb::Coord(j, k, 100))[i];
        }
    });
    saveData(data.constAccessor(), "src_#grid2.npy");
    
    SemiLagrangian3 solver;
    solver.advect(src, flow, 0.1, &dst, boundarySdf);
    
    data.forEachIndex([&](uint i, uint j, uint k) {
        if (i < 2) {
            data(i, j, k) = dst.getGrid()->tree().getValue(openvdb::Coord(j, k, 100))[i];
        }
    });
    saveData(data.constAccessor(), "dst_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SemiLagrangian3, Zalesak) {
    vox::Box3 box(vox::Vector3D(0.5 - 0.025, 0.6, 0.5),
                  vox::Vector3D(0.5 + 0.025, 0.85, 0.5));
    CellCenteredScalarGrid3 sdf(200, 200, 5,
                                1.0/200.0, 1.0/200.0, 1.0/5.0);
    CellCenteredScalarGrid3 sdf2(200, 200, 5,
                                1.0/200.0, 1.0/200.0, 1.0/5.0);
    sdf.getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                openvdb::Coord(199, 199, 4)),
                             0.0, true);
    sdf2.getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                 openvdb::Coord(199, 199, 4)),
                             0.0, true);
    
    sdf.fill([box](const vox::Vector3D& pt) {
        double disk = pt.distanceTo(vox::Vector3D(0.5, 0.75, 0.5)) - 0.15;
        double slot = box.closestDistance(pt);
        if (!box.boundingBox().contains(pt)) {
            slot *= -1.0;
        }
        return std::max(disk, slot);
    }, vox::ExecutionPolicy::kSerial);
    
    vox::CustomVectorField3 flow([](const vox::Vector3D& pt) {
        return vox::Vector3D(vox::kPiD / 3.14 * (0.5 - pt.y),
                             vox::kPiD / 3.14 * (pt.x - 0.5),
                             0.0);
    });
    
    vox::Array2<double> data(sdf.resolution().x, sdf.resolution().y);
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = sdf(openvdb::Coord(i, j, 2));
    });
    saveData(data.constAccessor(), "orig_#grid2,iso.npy");
    
    SemiLagrangian3 solver;
    
    for (int i = 0; i < 628; ++i) {
        solver.advect(sdf, flow, 0.02, &sdf2);
        sdf.swap(&sdf2);
    }
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = sdf(openvdb::Coord(i, j, 2));
    });
    saveData(data.constAccessor(), "rev0628_#grid2,iso.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SemiLagrangian3, Zalesak_cubic) {
    vox::Box3 box(vox::Vector3D(0.5 - 0.025, 0.6, 0.5),
                  vox::Vector3D(0.5 + 0.025, 0.85, 0.5));
    CellCenteredScalarGrid3 sdf(200, 200, 5,
                                1.0/200.0, 1.0/200.0, 1.0/5.0);
    CellCenteredScalarGrid3 sdf2(200, 200, 5,
                                1.0/200.0, 1.0/200.0, 1.0/5.0);
    sdf.getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                openvdb::Coord(199, 199, 4)),
                             0.0, true);
    sdf2.getGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                 openvdb::Coord(199, 199, 4)),
                             0.0, true);
    
    sdf.fill([box](const vox::Vector3D& pt) {
        double disk = pt.distanceTo(vox::Vector3D(0.5, 0.75, 0.5)) - 0.15;
        double slot = box.closestDistance(pt);
        if (!box.boundingBox().contains(pt)) {
            slot *= -1.0;
        }
        return std::max(disk, slot);
    }, vox::ExecutionPolicy::kSerial);
    
    vox::CustomVectorField3 flow([](const vox::Vector3D& pt) {
        return vox::Vector3D(vox::kPiD / 3.14 * (0.5 - pt.y),
                             vox::kPiD / 3.14 * (pt.x - 0.5),
                             0.0);
    });
    
    vox::Array2<double> data(sdf.resolution().x, sdf.resolution().y);
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = sdf(openvdb::Coord(i, j, 2));
    });
    saveData(data.constAccessor(), "orig_#grid2,iso.npy");
    
    CubicSemiLagrangian3 solver;
    
    for (int i = 0; i < 628; ++i) {
        solver.advect(sdf, flow, 0.02, &sdf2);
        sdf.swap(&sdf2);
    }
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = sdf(openvdb::Coord(i, j, 2));
    });
    saveData(data.constAccessor(), "rev0628_#grid2,iso.npy");
}
JET_END_TEST_F
