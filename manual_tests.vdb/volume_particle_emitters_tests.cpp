//
//  volume_particle_emitters_tests.cpp
//  manual_test.vdb
//
//  Created by Feng Yang on 2020/5/20.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "manual_tests.h"

#include "../src.common/sphere3.h"
#include "../src.common/surface_to_implicit3.h"
#include "../src.vdb/vdb_volume_particle_emitter3.hpp"
#include "../src.vdb/vdb_particle_system_solver3.hpp"

using namespace vdb;

JET_TESTS(VolumeParticleEmitter3);

JET_BEGIN_TEST_F(VolumeParticleEmitter3, EmitContinuousNonOverlapping) {
    ParticleSystemSolver3 solver;
    
    ParticleSystemData3Ptr particles = solver.particleSystemData();
    auto emitter
    = std::make_shared<VolumeParticleEmitter3>(
                                               std::make_shared<vox::SurfaceToImplicit3>(
                                                                                         std::make_shared<vox::Sphere3>(vox::Vector3D(), 1.0)),
                                               vox::BoundingBox3D(vox::Vector3D(-1, -1, -1), vox::Vector3D(1, 1, 1)),
                                               0.2);
    emitter->setIsOneShot(false);
    emitter->setAllowOverlapping(false);
    solver.setEmitter(emitter);
    
    saveParticleDataXy(particles, 0);
    
    for (vox::Frame frame(0, 1.0 / 60.0); frame.index < 120; ++frame) {
        solver.update(frame);
        
        saveParticleDataXy(particles, frame.index);
    }
}
JET_END_TEST_F
