//
//  particle_system_solver3_tests.cpp
//  manual_test.vdb
//
//  Created by Feng Yang on 2020/5/20.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "manual_tests.h"

#include "../src.common/rigid_body_collider3.h"
#include "../src.common/constant_vector_field3.h"
#include "../src.vdb/vdb_particle_system_solver3.hpp"
#include "../src.common/plane3.h"
#include "../src.vdb/vdb_point_particle_emitter3.hpp"

using namespace vdb;

JET_TESTS(ParticleSystemSolver3);

JET_BEGIN_TEST_F(ParticleSystemSolver3, PerfectBounce) {
    vox::Plane3Ptr plane = std::make_shared<vox::Plane3>(vox::Vector3D(0, 1, 0),
                                                         vox::Vector3D());
    vox::RigidBodyCollider3Ptr collider
    = std::make_shared<vox::RigidBodyCollider3>(plane);
    
    ParticleSystemSolver3 solver;
    solver.setCollider(collider);
    solver.setDragCoefficient(0.0);
    solver.setRestitutionCoefficient(1.0);
    
    ParticleSystemData3Ptr particles = solver.particleSystemData();
    particles->addParticle({0.0, 3.0, 0.0}, {1.0, 0.0, 0.0});
    
    vox::Array1<double> x(1000);
    vox::Array1<double> y(1000);
    char filename[256];
    snprintf(filename, sizeof(filename), "data.#line2,0000,x.npy");
    saveData(x.constAccessor(), 0, filename);
    snprintf(filename, sizeof(filename), "data.#line2,0000,y.npy");
    saveData(y.constAccessor(), 0, filename);
    
    vox::Frame frame;
    frame.timeIntervalInSeconds = 1.0 / 300.0;
    for (; frame.index < 1000; frame.advance()) {
        solver.update(frame);
        
        x[frame.index] = particles->positions()[0].x;
        y[frame.index] = particles->positions()[0].y;
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,x.npy",
                 frame.index);
        saveData(x.constAccessor(), frame.index, filename);
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,y.npy",
                 frame.index);
        saveData(y.constAccessor(), frame.index, filename);
    }
}
JET_END_TEST_F

JET_BEGIN_TEST_F(ParticleSystemSolver3, HalfBounce) {
    vox::Plane3Ptr plane = std::make_shared<vox::Plane3>(vox::Vector3D(0, 1, 0),
                                                         vox::Vector3D());
    vox::RigidBodyCollider3Ptr collider
    = std::make_shared<vox::RigidBodyCollider3>(plane);
    
    ParticleSystemSolver3 solver;
    solver.setCollider(collider);
    solver.setDragCoefficient(0.0);
    solver.setRestitutionCoefficient(0.5);
    
    ParticleSystemData3Ptr particles = solver.particleSystemData();
    particles->addParticle({0.0, 3.0, 0.0}, {1.0, 0.0, 0.0});
    
    vox::Array1<double> x(1000);
    vox::Array1<double> y(1000);
    char filename[256];
    snprintf(filename, sizeof(filename), "data.#line2,0000,x.npy");
    saveData(x.constAccessor(), 0, filename);
    snprintf(filename, sizeof(filename), "data.#line2,0000,y.npy");
    saveData(y.constAccessor(), 0, filename);
    
    vox::Frame frame;
    frame.timeIntervalInSeconds = 1.0 / 300.0;
    for (; frame.index < 1000; frame.advance()) {
        solver.update(frame);
        
        x[frame.index] = particles->positions()[0].x;
        y[frame.index] = particles->positions()[0].y;
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,x.npy",
                 frame.index);
        saveData(x.constAccessor(), frame.index, filename);
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,y.npy",
                 frame.index);
        saveData(y.constAccessor(), frame.index, filename);
    }
}
JET_END_TEST_F

JET_BEGIN_TEST_F(ParticleSystemSolver3, HalfBounceWithFriction) {
    vox::Plane3Ptr plane = std::make_shared<vox::Plane3>(vox::Vector3D(0, 1, 0),
                                                         vox::Vector3D());
    vox::RigidBodyCollider3Ptr collider
    = std::make_shared<vox::RigidBodyCollider3>(plane);
    collider->setFrictionCoefficient(0.04);
    
    ParticleSystemSolver3 solver;
    solver.setCollider(collider);
    solver.setDragCoefficient(0.0);
    solver.setRestitutionCoefficient(0.5);
    
    ParticleSystemData3Ptr particles = solver.particleSystemData();
    particles->addParticle({0.0, 3.0, 0.0}, {1.0, 0.0, 0.0});
    
    vox::Array1<double> x(1000);
    vox::Array1<double> y(1000);
    char filename[256];
    snprintf(filename, sizeof(filename), "data.#line2,0000,x.npy");
    saveData(x.constAccessor(), 0, filename);
    snprintf(filename, sizeof(filename), "data.#line2,0000,y.npy");
    saveData(y.constAccessor(), 0, filename);
    
    vox::Frame frame;
    frame.timeIntervalInSeconds = 1.0 / 300.0;
    for (; frame.index < 1000; frame.advance()) {
        solver.update(frame);
        
        x[frame.index] = particles->positions()[0].x;
        y[frame.index] = particles->positions()[0].y;
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,x.npy",
                 frame.index);
        saveData(x.constAccessor(), frame.index, filename);
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,y.npy",
                 frame.index);
        saveData(y.constAccessor(), frame.index, filename);
    }
}
JET_END_TEST_F

JET_BEGIN_TEST_F(ParticleSystemSolver3, NoBounce) {
    vox::Plane3Ptr plane = std::make_shared<vox::Plane3>(vox::Vector3D(0, 1, 0),
                                                         vox::Vector3D());
    vox::RigidBodyCollider3Ptr collider
    = std::make_shared<vox::RigidBodyCollider3>(plane);
    
    ParticleSystemSolver3 solver;
    solver.setCollider(collider);
    solver.setDragCoefficient(0.0);
    solver.setRestitutionCoefficient(0.0);
    
    ParticleSystemData3Ptr particles = solver.particleSystemData();
    particles->addParticle({0.0, 3.0, 0.0}, {1.0, 0.0, 0.0});
    
    vox::Array1<double> x(1000);
    vox::Array1<double> y(1000);
    char filename[256];
    snprintf(filename, sizeof(filename), "data.#line2,0000,x.npy");
    saveData(x.constAccessor(), 0, filename);
    snprintf(filename, sizeof(filename), "data.#line2,0000,y.npy");
    saveData(y.constAccessor(), 0, filename);
    
    vox::Frame frame;
    frame.timeIntervalInSeconds = 1.0 / 300.0;
    for (; frame.index < 1000; frame.advance()) {
        solver.update(frame);
        
        x[frame.index] = particles->positions()[0].x;
        y[frame.index] = particles->positions()[0].y;
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,x.npy",
                 frame.index);
        saveData(x.constAccessor(), frame.index, filename);
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#line2,%04d,y.npy",
                 frame.index);
        saveData(y.constAccessor(), frame.index, filename);
    }
}
JET_END_TEST_F

JET_BEGIN_TEST_F(ParticleSystemSolver3, Update) {
    auto plane = vox::Plane3::builder()
    .withNormal({0, 1, 0})
    .withPoint({0, 0, 0})
    .makeShared();
    
    auto collider = vox::RigidBodyCollider3::builder()
    .withSurface(plane)
    .makeShared();
    
    auto wind = vox::ConstantVectorField3::builder()
    .withValue({1, 0, 0})
    .makeShared();
    
    auto emitter = PointParticleEmitter3::builder()
    .withOrigin({0, 3, 0})
    .withDirection({0, 1, 0})
    .withSpeed(5)
    .withSpreadAngleInDegrees(45.0)
    .withMaxNumberOfNewParticlesPerSecond(300)
    .makeShared();
    
    auto solver = ParticleSystemSolver3::builder().makeShared();
    solver->setCollider(collider);
    solver->setEmitter(emitter);
    solver->setWind(wind);
    solver->setDragCoefficient(0.0);
    solver->setRestitutionCoefficient(0.5);
    
    for (vox::Frame frame(0, 1.0 / 60.0); frame.index < 360; ++frame) {
        solver->update(frame);
        
        saveParticleDataXy(solver->particleSystemData(), frame.index);
    }
}
JET_END_TEST_F
