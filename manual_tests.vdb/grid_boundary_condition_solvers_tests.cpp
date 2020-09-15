//
//  grid_boundary_condition_solvers_tests.cpp
//  vdb_manual_tests
//
//  Created by Feng Yang on 2020/2/12.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "manual_tests.h"

#include "../src.common/array2.h"
#include "../src.vdb/vdb_blocked_boundary_condition_solver3.hpp"
#include "../src.vdb/vdb_fractional_boundary_condition_solver3.hpp"
#include "../src.common/implicit_surface_set3.h"
#include "../src.common/plane3.h"
#include "../src.common/rigid_body_collider3.h"
#include "../src.common/sphere3.h"
#include "../src.common/surface_to_implicit3.h"
#include "../src.common/triangle_mesh3.h"

using namespace vdb;

JET_TESTS(GridBlockedBoundaryConditionSolver3);

JET_BEGIN_TEST_F(GridBlockedBoundaryConditionSolver3, ConstrainVelocitySmall) {
    BlockedBoundaryConditionSolver3 solver;
    auto collider
    = std::make_shared<vox::RigidBodyCollider3>(
                                                std::make_shared<vox::Plane3>(vox::Vector3D(1, 1, 0).normalized(),
                                                                              vox::Vector3D()));
    vox::Size3 gridSize(10, 10, 10);
    vox::Vector3D gridSpacing(1.0, 1.0, 1.0);
    vox::Vector3D gridOrigin(-5.0, -5.0, -5.0);
    
    solver.updateCollider(collider, gridSize, gridSpacing, gridOrigin);
    FaceCenteredGrid3 velocity(gridSize, gridSpacing, gridOrigin);
    velocity.getUGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(10, 9, 9)), 0.0, true);
    velocity.getVGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(9, 10, 9)), 0.0, true);
    velocity.getWGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(9, 9, 10)), 0.0, true);
    velocity.fill(vox::Vector3D(1.0, 1.0, 0.0), vox::ExecutionPolicy::kSerial);
    
    solver.constrainVelocity(&velocity);
    
    // Output
    vox::Array2<double> dataU(10, 10);
    vox::Array2<double> dataV(10, 10);
    vox::Array2<double> dataM(10, 10);
    
    dataU.forEachIndex([&](uint i, uint j) {
        vox::Vector3D vel = velocity.valueAtCellCenter(i, j, 5);
        dataU(i, j) = vel.x;
        dataV(i, j) = vel.y;
        dataM(i, j) = static_cast<double>(solver.marker()(i, j, 5));
    });
    
    saveData(dataU.constAccessor(), "data_#grid2,x.npy");
    saveData(dataV.constAccessor(), "data_#grid2,y.npy");
    saveData(dataM.constAccessor(), "marker_#grid2.npy");
    
}
JET_END_TEST_F

JET_BEGIN_TEST_F(GridBlockedBoundaryConditionSolver3, ConstrainVelocity) {
    double dx = 1.0 / 32.0;
    FaceCenteredGrid3 velocity(vox::Size3(64, 32, 32),
                               vox::Vector3D(dx, dx, dx),
                               vox::Vector3D());
    velocity.getUGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(64, 31, 31)), 0.0, true);
    velocity.getVGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(63, 32, 31)), 0.0, true);
    velocity.getWGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(63, 31, 32)), 0.0, true);
    velocity.fill(vox::Vector3D(1.0, 0.0, 0.0), vox::ExecutionPolicy::kSerial);
    
    // Collider setting
    auto surfaceSet = std::make_shared<vox::ImplicitSurfaceSet3>();
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(1, 0.5, 0.5), 0.25));
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(2,1,0.5), 0.25));
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(0,0,0.5), 0.25));
    vox::RigidBodyCollider3Ptr collider
    = std::make_shared<vox::RigidBodyCollider3>(surfaceSet);
    collider->linearVelocity = vox::Vector3D(-1.0, 0.0, 0.0);
    
    // Solver setting
    BlockedBoundaryConditionSolver3 solver;
    solver.updateCollider(collider,
                          velocity.resolution(),
                          velocity.gridSpacing(),
                          velocity.origin());
    solver.setClosedDomainBoundaryFlag(
                                       vox::kDirectionRight |
                                       vox::kDirectionDown |
                                       vox::kDirectionUp);
    
    // Constrain velocity
    solver.constrainVelocity(&velocity, 5);
    
    // Output
    vox::Array2<double> dataU(64, 32);
    vox::Array2<double> dataV(64, 32);
    vox::Array2<double> dataM(64, 32);
    dataU.forEachIndex([&](uint i, uint j) {
        vox::Vector3D vel = velocity.valueAtCellCenter(i, j, 16);
        dataU(i, j) = vel.x;
        dataV(i, j) = vel.y;
        dataM(i, j) = static_cast<double>(solver.marker()(i, j, 16));
    });
    
    saveData(dataU.constAccessor(), "data_#grid2,x.npy");
    saveData(dataV.constAccessor(), "data_#grid2,y.npy");
    saveData(dataM.constAccessor(), "marker_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(
                 GridBlockedBoundaryConditionSolver3, ConstrainVelocityWithFriction) {
    double dx = 1.0 / 32.0;
    FaceCenteredGrid3 velocity(vox::Size3(64, 32, 32),
                               vox::Vector3D(dx, dx, dx),
                               vox::Vector3D());
    velocity.getUGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(64, 31, 31)), 0.0, true);
    velocity.getVGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(63, 32, 31)), 0.0, true);
    velocity.getWGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(63, 31, 32)), 0.0, true);
    velocity.fill(vox::Vector3D(1.0, 0.0, 0.0), vox::ExecutionPolicy::kSerial);
    
    // Collider setting
    auto surfaceSet = std::make_shared<vox::ImplicitSurfaceSet3>();
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(1, 0.5, 0.5), 0.25));
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(2,1,0.5), 0.25));
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(0,0,0.5), 0.25));
    vox::RigidBodyCollider3Ptr collider
    = std::make_shared<vox::RigidBodyCollider3>(surfaceSet);
    collider->linearVelocity = vox::Vector3D(-1.0, 0.0, 0.0);
    collider->setFrictionCoefficient(1.0);
    
    // Solver setting
    BlockedBoundaryConditionSolver3 solver;
    solver.updateCollider(
                          collider,
                          velocity.resolution(),
                          velocity.gridSpacing(),
                          velocity.origin());
    solver.setClosedDomainBoundaryFlag(
                                       vox::kDirectionRight |
                                       vox::kDirectionDown |
                                       vox::kDirectionUp);
    
    // Constrain velocity
    solver.constrainVelocity(&velocity, 5);
    
    // Output
    vox::Array2<double> dataU(64, 32);
    vox::Array2<double> dataV(64, 32);
    vox::Array2<double> dataM(64, 32);
    dataU.forEachIndex([&](uint i, uint j) {
        vox::Vector3D vel = velocity.valueAtCellCenter(i, j, 16);
        dataU(i, j) = vel.x;
        dataV(i, j) = vel.y;
        dataM(i, j) = static_cast<double>(solver.marker()(i, j, 16));
    });
    
    saveData(dataU.constAccessor(), "data_#grid2,x.npy");
    saveData(dataV.constAccessor(), "data_#grid2,y.npy");
    saveData(dataM.constAccessor(), "marker_#grid2.npy");
}
JET_END_TEST_F

//-------------------------------------------------------------------------
JET_TESTS(GridFractionalBoundaryConditionSolver3);

JET_BEGIN_TEST_F(
                 GridFractionalBoundaryConditionSolver3,
                 ConstrainVelocity) {
    double dx = 1.0 / 32.0;
    FaceCenteredGrid3 velocity(vox::Size3(64, 32, 32),
                               vox::Vector3D(dx, dx, dx),
                               vox::Vector3D());
    velocity.getUGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(64, 31, 31)), 0.0, true);
    velocity.getVGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(63, 32, 31)), 0.0, true);
    velocity.getWGrid()->denseFill(openvdb::CoordBBox(openvdb::Coord(0, 0, 0),
                                                      openvdb::Coord(63, 31, 32)), 0.0, true);
    velocity.fill(vox::Vector3D(1.0, 0.0, 0.0), vox::ExecutionPolicy::kSerial);
    
    // Collider setting
    auto surfaceSet = std::make_shared<vox::ImplicitSurfaceSet3>();
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(1, 0.5, 0.5), 0.25));
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(2,1,0.5), 0.25));
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(vox::Vector3D(0,0,0.5), 0.25));
    vox::RigidBodyCollider3Ptr collider
    = std::make_shared<vox::RigidBodyCollider3>(surfaceSet);
    collider->linearVelocity = vox::Vector3D(-1.0, 0.0, 0.0);
    
    // Solver setting
    FractionalBoundaryConditionSolver3 solver;
    solver.updateCollider(collider,
                          velocity.resolution(),
                          velocity.gridSpacing(),
                          velocity.origin());
    solver.setClosedDomainBoundaryFlag(
                                       vox::kDirectionRight |
                                       vox::kDirectionDown |
                                       vox::kDirectionUp);
    
    // Constrain velocity
    solver.constrainVelocity(&velocity, 5);
    
    // Output
    vox::Array2<double> dataU(64, 32);
    vox::Array2<double> dataV(64, 32);
    
    dataU.forEachIndex([&](uint i, uint j) {
        vox::Vector3D vel = velocity.valueAtCellCenter(i, j, 16);
        dataU(i, j) = vel.x;
        dataV(i, j) = vel.y;
    });
    
    saveData(dataU.constAccessor(), "data_#grid2,x.npy");
    saveData(dataV.constAccessor(), "data_#grid2,y.npy");
}
JET_END_TEST_F
