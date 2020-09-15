// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "manual_tests.h"

#include "../src.common/array3.h"
#include "../src.common/marching_cubes.h"
#include "../src.vdb/vdb_spherical_points_to_implicit3.h"
#include "../src.vdb/vdb_vertex_centered_scalar_grid3.h"
#include "../src.common/triangle_mesh3.h"
#include <openvdb/tools/VolumeToMesh.h>

#include <random>

using namespace vdb;

JET_TESTS(SphericalPointsToImplicit3);

JET_BEGIN_TEST_F(SphericalPointsToImplicit3, ConvertTwo) {
    vox::Array1<vox::Vector3D> points;
    
    std::mt19937 rng{0};
    std::uniform_real_distribution<> dist(0.2, 0.8);
    for (size_t i = 0; i < 2; ++i) {
        points.append({dist(rng), dist(rng), dist(rng)});
    }
    
    VertexCenteredScalarGrid3 grid(128, 128, 128,
                                   1.0 / 128, 1.0 / 128, 1.0 / 128);
    
    SphericalPointsToImplicit3 converter(0.15);
    converter.convert(points, &grid);
    
    vox::Array3<double> val(grid.dataSize());
    grid.getData(val.accessor());
    
    vox::TriangleMesh3 triMesh;
    marchingCubes(val,
                  grid.gridSpacing(),
                  grid.dataOrigin(), &triMesh, 0, vox::kDirectionAll);
    
    saveTriangleMeshData(triMesh,
                         "spherical_points_to_implicit3_convert_two.obj");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphericalPointsToImplicit3, ConvertMany) {
    vox::Array1<vox::Vector3D> points;
    
    std::mt19937 rng{0};
    std::uniform_real_distribution<> dist(0.2, 0.8);
    for (size_t i = 0; i < 500; ++i) {
        points.append({dist(rng), dist(rng), dist(rng)});
    }
    
    VertexCenteredScalarGrid3 grid(128, 128, 128,
                                   1.0 / 128, 1.0 / 128, 1.0 / 128);
    
    SphericalPointsToImplicit3 converter(0.1);
    converter.convert(points, &grid);
    
    vox::Array3<double> val(grid.dataSize());
    grid.getData(val.accessor());
    
    vox::TriangleMesh3 triMesh;
    marchingCubes(val,
                  grid.gridSpacing(),
                  grid.dataOrigin(), &triMesh, 0, vox::kDirectionAll);
    
    saveTriangleMeshData(triMesh,
                         "spherical_points_to_implicit3_convert_many.obj");
}
JET_END_TEST_F
