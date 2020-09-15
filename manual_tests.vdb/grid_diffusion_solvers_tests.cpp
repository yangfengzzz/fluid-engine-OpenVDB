// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "manual_tests.h"

#include "../src.common/array2.h"
#include "../src.vdb/vdb_backward_euler_diffusion_solver3.hpp"
#include "../src.vdb/vdb_forward_euler_diffusion_solver3.hpp"
#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.common/custom_scalar_field3.h"


using namespace vdb;

JET_TESTS(GridForwardEulerDiffusionSolver3);

JET_BEGIN_TEST_F(GridForwardEulerDiffusionSolver3, Solve) {
    vox::Size3 size(160, 120, 150);
    vox::Vector3D gridSpacing(1.0/size.x, 1.0/size.x, 1.0/size.x);
    
    CellCenteredScalarGrid3 src(size, gridSpacing);
    CellCenteredScalarGrid3 dst(size, gridSpacing);
    vox::Array2<double> data(160, 120);
    
    src.fill([&](const vox::Vector3D& x)->double {
        return ((x-src.boundingBox().midPoint()).length() < 0.2) ? 1.0 : 0.0;
    }, vox::ExecutionPolicy::kSerial);
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src(openvdb::Coord(i, j, 75));
    });
    
    double timeStep = 0.01;
    double diffusionCoeff = vox::square(gridSpacing.x) / timeStep / 12.0;
    
    ForwardEulerDiffusionSolver3 diffusionSolver;
    
    diffusionSolver.solve(src, diffusionCoeff, timeStep, &dst);
    dst.swap(&src);
    
    saveData(data.constAccessor(), "src_#grid2.npy");
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src.getGrid()->tree().getValue(openvdb::Coord(i, j, 75) );
    });
    
    saveData(data.constAccessor(), "dst_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(GridForwardEulerDiffusionSolver3, Unstable) {
    vox::Size3 size(160, 120, 150);
    vox::Vector3D gridSpacing(1.0/size.x, 1.0/size.x, 1.0/size.x);
    
    CellCenteredScalarGrid3 src(size, gridSpacing);
    CellCenteredScalarGrid3 dst(size, gridSpacing);
    src.fill([&](const vox::Vector3D& x) {
        return ((x- src.boundingBox().midPoint()).length() < 0.2) ? 1.0 : 0.0;
    }, vox::ExecutionPolicy::kSerial);
    
    vox::Array2<double> data(160, 120);
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src(openvdb::Coord(i, j, 75) );
    });
    
    double timeStep = 0.01;
    double diffusionCoeff = vox::square(gridSpacing.x) / timeStep / 12.0;
    
    ForwardEulerDiffusionSolver3 diffusionSolver;
    
    diffusionSolver.solve(
                          src,
                          10.0 * diffusionCoeff,
                          timeStep,
                          &dst);
    dst.swap(&src);
    
    saveData(data.constAccessor(), "src_#grid2.npy");
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src.getGrid()->tree().getValue(openvdb::Coord(i, j, 75) );
    });
    
    saveData(data.constAccessor(), "dst_#grid2.npy");
}
JET_END_TEST_F


JET_TESTS(GridBackwardEulerDiffusionSolver3);

JET_BEGIN_TEST_F(GridBackwardEulerDiffusionSolver3, Solve) {
    vox::Size3 size(160, 120, 150);
    vox::Vector3D gridSpacing(1.0/size.x, 1.0/size.x, 1.0/size.x);
    
    CellCenteredScalarGrid3 src(size, gridSpacing);
    CellCenteredScalarGrid3 dst(size, gridSpacing);
    src.fill([&](const vox::Vector3D& x) {
        return ((x- src.boundingBox().midPoint()).length() < 0.2) ? 1.0 : 0.0;
    }, vox::ExecutionPolicy::kSerial);
    
    vox::Array2<double> data(160, 120);
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src.getGrid()->tree().getValue(openvdb::Coord(i, j, 75) );
    });
    
    double timeStep = 0.01;
    double diffusionCoeff = vox::square(gridSpacing.x) / timeStep / 12.0;
    
    BackwardEulerDiffusionSolver3 diffusionSolver;
    
    diffusionSolver.solve(src, diffusionCoeff, timeStep, &dst);
    dst.swap(&src);
    
    saveData(data.constAccessor(), "src_#grid2.npy");
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src.getGrid()->tree().getValue(openvdb::Coord(i, j, 75) );
    });
    
    saveData(data.constAccessor(), "dst_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(GridBackwardEulerDiffusionSolver3, Stable) {
    vox::Size3 size(160, 120, 150);
    vox::Vector3D gridSpacing(1.0/size.x, 1.0/size.x, 1.0/size.x);
    
    CellCenteredScalarGrid3 src(size, gridSpacing);
    CellCenteredScalarGrid3 dst(size, gridSpacing);
    src.fill([&](const vox::Vector3D& x) {
        return ((x- src.boundingBox().midPoint()).length() < 0.2) ? 1.0 : 0.0;
    }, vox::ExecutionPolicy::kSerial);
    
    vox::Array2<double> data(160, 120);
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src.getGrid()->tree().getValue(openvdb::Coord(i, j, 75) );
    });
    
    double timeStep = 0.01;
    double diffusionCoeff = vox::square(gridSpacing.x) / timeStep / 12.0;
    
    BackwardEulerDiffusionSolver3 diffusionSolver;
    
    diffusionSolver.solve(
                          src,
                          10.0 * diffusionCoeff,
                          timeStep,
                          &dst);
    dst.swap(&src);
    
    saveData(data.constAccessor(), "src_#grid2.npy");
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src.getGrid()->tree().getValue(openvdb::Coord(i, j, 75) );
    });
    
    saveData(data.constAccessor(), "dst_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(
                 GridBackwardEulerDiffusionSolver3,
                 SolveWithBoundaryDirichlet) {
    vox::Size3 size(80, 60, 75);
    vox::Vector3D gridSpacing(1.0/size.x, 1.0/size.x, 1.0/size.x);
    
    CellCenteredScalarGrid3 src(size, gridSpacing);
    CellCenteredScalarGrid3 dst(size, gridSpacing);
    
    vox::Vector3D boundaryCenter = src.boundingBox().midPoint();
    vox::CustomScalarField3 boundarySdf(
                                        [&](const vox::Vector3D& x) {
        return boundaryCenter.x - x.x;
    });
    
    vox::Array2<double> data(size.x, size.y);
    
    
    src.fill([&](const vox::Vector3D& x) {
        return ((x- src.boundingBox().midPoint()).length() < 0.2) ? 1.0 : 0.0;
    }, vox::ExecutionPolicy::kSerial);
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src(openvdb::Coord(i, j, (uint)size.z / 2) );
    });
    
    double timeStep = 0.01;
    double diffusionCoeff = 100 * vox::square(gridSpacing.x) / timeStep / 12.0;
    
    BackwardEulerDiffusionSolver3 diffusionSolver(
                                                  BackwardEulerDiffusionSolver3::Dirichlet);
    
    diffusionSolver.solve(src, diffusionCoeff, timeStep, &dst, boundarySdf);
    dst.swap(&src);
    
    saveData(data.constAccessor(), "src_#grid2.npy");
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src(openvdb::Coord(i, j, (uint)size.z / 2) );
    });
    
    saveData(data.constAccessor(), "dst_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(GridBackwardEulerDiffusionSolver3, SolveWithBoundaryNeumann) {
    vox::Size3 size(80, 60, 75);
    vox::Vector3D gridSpacing(1.0/size.x, 1.0/size.x, 1.0/size.x);
    
    CellCenteredScalarGrid3 src(size, gridSpacing);
    CellCenteredScalarGrid3 dst(size, gridSpacing);
    
    vox::Vector3D boundaryCenter = src.boundingBox().midPoint();
    vox::CustomScalarField3 boundarySdf(
                                        [&](const vox::Vector3D& x) {
        return boundaryCenter.x - x.x;
    });
    
    vox::Array2<double> data(size.x, size.y);
    
    src.fill([&](const vox::Vector3D& x) {
        return ((x- src.boundingBox().midPoint()).length() < 0.2) ? 1.0 : 0.0;
    }, vox::ExecutionPolicy::kSerial);
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src(openvdb::Coord(i, j, (uint)size.z / 2) );
    });
    
    double timeStep = 0.01;
    double diffusionCoeff = 100 * vox::square(gridSpacing.x) / timeStep / 12.0;
    
    BackwardEulerDiffusionSolver3 diffusionSolver(
                                                  BackwardEulerDiffusionSolver3::Neumann);
    
    diffusionSolver.solve(src, diffusionCoeff, timeStep, &dst, boundarySdf);
    dst.swap(&src);
    
    saveData(data.constAccessor(), "src_#grid2.npy");
    
    data.forEachIndex([&](uint i, uint j) {
        data(i, j) = src(openvdb::Coord(i, j, (uint)size.z / 2) );
    });
    
    saveData(data.constAccessor(), "dst_#grid2.npy");
}
JET_END_TEST_F
