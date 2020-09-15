//
//  sph_system_data_tests.cpp
//  manual_test.vdb
//
//  Created by Feng Yang on 2020/5/20.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "manual_tests.h"

#include "../src.common/bcc_lattice_point_generator.h"
#include "../src.common/bounding_box3.h"
#include "../src.common/cell_centered_scalar_grid2.h"
#include "../src.common/triangle_point_generator.h"
#include "../src.common/parallel.h"
#include "../src.vdb/vdb_sph_system_data3.hpp"
#include <random>

using namespace vdb;

JET_TESTS(SphSystemData3);

JET_BEGIN_TEST_F(SphSystemData3, Interpolate) {
    vox::Array1<vox::Vector3D> points;
    vox::BccLatticePointGenerator pointsGenerator;
    vox::BoundingBox3D bbox(
                            vox::Vector3D(0, 0, 0),
                            vox::Vector3D(1, 1, 1));
    double spacing = 0.1;
    
    pointsGenerator.generate(bbox, spacing, &points);
    
    SphSystemData3 sphSystem;
    sphSystem.addParticles(
                           vox::ConstArrayAccessor1<vox::Vector3D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();
    
    vox::Array1<double> data(points.size(), 1.0);
    
    vox::CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);
    
    auto gridPos = grid.dataPosition();
    vox::parallelFor(vox::kZeroSize, grid.dataSize().x,
                     vox::kZeroSize, grid.dataSize().y,
                     [&](size_t i, size_t j) {
        vox::Vector2D xy = gridPos(i, j);
        vox::Vector3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, data);
    });
    
    saveData(grid.constDataAccessor(), "data_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphSystemData3, Gradient) {
    vox::Array1<vox::Vector3D> points;
    vox::BccLatticePointGenerator pointsGenerator;
    vox::BoundingBox3D bbox(
                            vox::Vector3D(0, 0, 0),
                            vox::Vector3D(1, 1, 1));
    double spacing = 0.1;
    
    pointsGenerator.generate(bbox, spacing, &points);
    
    SphSystemData3 sphSystem;
    sphSystem.addParticles(
                           vox::ConstArrayAccessor1<vox::Vector3D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();
    
    vox::Array1<double> data(points.size());
    vox::Array1<double> gradX(points.size()), gradY(points.size());
    std::mt19937 rng(0);
    std::uniform_real_distribution<> d(0.0, 1.0);
    
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = d(rng);
    }
    
    for (size_t i = 0; i < data.size(); ++i) {
        vox::Vector3D g = sphSystem.gradientAt(i, data);
        gradX[i] = g.x;
        gradY[i] = g.y;
    }
    
    vox::CellCenteredScalarGrid2 grid(64, 64, 1.0 / 64, 1.0 / 64);
    vox::CellCenteredScalarGrid2 grid2(64, 64, 1.0 / 64, 1.0 / 64);
    
    auto gridPos = grid.dataPosition();
    vox::parallelFor(
                     vox::kZeroSize,
                     grid.dataSize().x,
                     vox::kZeroSize,
                     grid.dataSize().y,
                     [&](size_t i, size_t j) {
        vox::Vector2D xy = gridPos(i, j);
        vox::Vector3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, data);
    });
    
    saveData(grid.constDataAccessor(), "data_#grid2.npy");
    
    vox::parallelFor(vox::kZeroSize, grid.dataSize().x,
                     vox::kZeroSize, grid.dataSize().y,
                     [&](size_t i, size_t j) {
        vox::Vector2D xy = gridPos(i, j);
        vox::Vector3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, gradX);
        grid2(i, j) = sphSystem.interpolate(p, gradY);
    });
    
    saveData(grid.constDataAccessor(), "gradient_#grid2,x.npy");
    saveData(grid2.constDataAccessor(), "gradient_#grid2,y.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphSystemData3, Laplacian) {
    vox::Array1<vox::Vector3D> points;
    vox::BccLatticePointGenerator pointsGenerator;
    vox::BoundingBox3D bbox(
                            vox::Vector3D(0, 0, 0),
                            vox::Vector3D(1, 1, 1));
    double spacing = 0.1;
    
    pointsGenerator.generate(bbox, spacing, &points);
    
    SphSystemData3 sphSystem;
    sphSystem.addParticles(
                           vox::ConstArrayAccessor1<vox::Vector3D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();
    
    vox::Array1<double> data(points.size()), laplacian(points.size());
    std::mt19937 rng(0);
    std::uniform_real_distribution<> d(0.0, 1.0);
    
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = d(rng);
    }
    
    for (size_t i = 0; i < data.size(); ++i) {
        laplacian[i] = sphSystem.laplacianAt(i, data);
    }
    
    vox::CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);
    
    auto gridPos = grid.dataPosition();
    vox::parallelFor(vox::kZeroSize, grid.dataSize().x,
                     vox::kZeroSize, grid.dataSize().y,
                     [&](size_t i, size_t j) {
        vox::Vector2D xy = gridPos(i, j);
        vox::Vector3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, data);
    });
    
    saveData(grid.constDataAccessor(), "data_#grid2.npy");
    
    vox::parallelFor(vox::kZeroSize, grid.dataSize().x,
                     vox::kZeroSize, grid.dataSize().y,
                     [&](size_t i, size_t j) {
        vox::Vector2D xy = gridPos(i, j);
        vox::Vector3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, laplacian);
    });
    
    saveData(grid.constDataAccessor(), "laplacian_#grid2.npy");
}
JET_END_TEST_F
