//
//  sph_solver3_tests.cpp
//  manual_test.vdb
//
//  Created by Feng Yang on 2020/5/20.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "manual_tests.h"

#include "../src.common/box3.h"
#include "../src.common/implicit_surface_set3.h"
#include "../src.common/plane3.h"
#include "../src.common/rigid_body_collider3.h"
#include "../src.vdb/vdb_sph_solver3.hpp"
#include "../src.common/sphere3.h"
#include "../src.common/surface_to_implicit3.h"
#include "../src.vdb/vdb_volume_particle_emitter3.hpp"

using namespace vdb;

JET_TESTS(SphSolver3);

JET_BEGIN_TEST_F(SphSolver3, SteadyState) {
    SphSolver3 solver;
    solver.setViscosityCoefficient(0.1);
    solver.setPseudoViscosityCoefficient(10.0);
    
    SphSystemData3Ptr particles = solver.sphSystemData();
    const double targetSpacing = particles->targetSpacing();
    
    vox::BoundingBox3D initialBound(vox::Vector3D(), vox::Vector3D(1, 0.5, 1));
    initialBound.expand(-targetSpacing);
    
    auto emitter
    = std::make_shared<VolumeParticleEmitter3>(
                                               std::make_shared<vox::SurfaceToImplicit3>(
                                                                                         std::make_shared<vox::Sphere3>(vox::Vector3D(), 10.0)),
                                               initialBound,
                                               particles->targetSpacing(),
                                               vox::Vector3D());
    emitter->setJitter(0.0);
    solver.setEmitter(emitter);
    
    vox::Box3Ptr box = std::make_shared<vox::Box3>(vox::Vector3D(), vox::Vector3D(1, 1, 1));
    box->isNormalFlipped = true;
    vox::RigidBodyCollider3Ptr collider = std::make_shared<vox::RigidBodyCollider3>(box);
    
    solver.setCollider(collider);
    
    saveParticleDataXy(particles, 0);
    
    for (vox::Frame frame(0, 1.0 / 60.0); frame.index < 100; frame.advance()) {
        solver.update(frame);
        
        saveParticleDataXy(particles, frame.index);
    }
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphSolver3, WaterDrop) {
    const double targetSpacing = 0.02;
    
    vox::BoundingBox3D domain(vox::Vector3D(), vox::Vector3D(1, 2, 0.5));
    
    // Initialize solvers
    SphSolver3 solver;
    solver.setPseudoViscosityCoefficient(0.0);
    
    SphSystemData3Ptr particles = solver.sphSystemData();
    particles->setTargetDensity(1000.0);
    particles->setTargetSpacing(targetSpacing);
    
    // Initialize source
    vox::ImplicitSurfaceSet3Ptr surfaceSet = std::make_shared<vox::ImplicitSurfaceSet3>();
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Plane3>(
                                                                 vox::Vector3D(0, 1, 0),
                                                                 vox::Vector3D(0, 0.25 * domain.height(), 0)));
    surfaceSet->addExplicitSurface(
                                   std::make_shared<vox::Sphere3>(
                                                                  domain.midPoint(), 0.15 * domain.width()));
    
    vox::BoundingBox3D sourceBound(domain);
    sourceBound.expand(-targetSpacing);
    
    auto emitter = std::make_shared<VolumeParticleEmitter3>(
                                                            surfaceSet,
                                                            sourceBound,
                                                            targetSpacing,
                                                            vox::Vector3D());
    solver.setEmitter(emitter);
    
    // Initialize boundary
    vox::Box3Ptr box = std::make_shared<vox::Box3>(domain);
    box->isNormalFlipped = true;
    vox::RigidBodyCollider3Ptr collider = std::make_shared<vox::RigidBodyCollider3>(box);
    solver.setCollider(collider);
    
    // Make it fast, but stable
    solver.setViscosityCoefficient(0.01);
    solver.setTimeStepLimitScale(5.0);
    
    saveParticleDataXy(particles, 0);
    
    for (vox::Frame frame(0, 1.0 / 60.0); frame.index < 100; frame.advance()) {
        solver.update(frame);
        
        saveParticleDataXy(particles, frame.index);
    }
}
JET_END_TEST_F

