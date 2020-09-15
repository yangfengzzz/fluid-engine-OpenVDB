// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "manual_tests.h"

#include "../src.common/array2.h"
#include "../src.common/box3.h"
#include "../src.vdb/vdb_fractional_boundary_condition_solver3.hpp"
#include "../src.vdb/vdb_fractional_single_phase_pressure_solver3.hpp"
#include "../src.vdb/vdb_smoke_solver3.hpp"
#include "../src.common/implicit_surface_set3.h"
#include "../src.common/level_set_utils.h"
#include "../src.common/rigid_body_collider3.h"
#include "../src.common/sphere3.h"
#include "../src.common/surface_to_implicit3.h"
#include "../src.vdb/vdb_volume_grid_emitter3.hpp"
#include "../src.vdb/vdb_volume_particle_emitter3.hpp"
#include "../src.common/triangle_mesh3.h"
#include <algorithm>

using namespace vdb;

/*
JET_TESTS(GridSmokeSolver3);

JET_BEGIN_TEST_F(GridSmokeSolver3, Rising) {
    //
    // This is a replica of smoke_sim example 4
    //
    
    uint resolutionX = 50;
    
    // Build solver
    auto solver = GridSmokeSolver3::builder()
    .withResolution({resolutionX, 6 * resolutionX / 5, resolutionX / 2})
    .withDomainSizeX(1.0)
    .makeShared();
    
    solver->setBuoyancyTemperatureFactor(2.0);
    
    // Build emitter
    auto box = vox::Box3::builder()
    .withLowerCorner({0.05, 0.1, 0.225})
    .withUpperCorner({0.1, 0.15, 0.275})
    .makeShared();
    
    auto particle_emitter = VolumeParticleEmitter3::builder()
    .withSurface(box)
    .withSpacing(1.0 / 1024.0)
    .withIsOneShot(true)
    .makeShared();
    solver->setParticleEmitter(particle_emitter);
    
    auto emitter = GridEmitter3::builder()
    .withSourceRegion(box)
    .withIsOneShot(false)
    .makeShared();
    
    solver->setEmitter(emitter);
    emitter->addStepFunctionTarget(solver->smokeDensity(), 0, 1);
    emitter->addStepFunctionTarget(solver->temperature(), 0, 1);
    emitter->addTarget(
                       solver->velocity(),
                       [](double sdf,
                          const vox::Vector3D& pt,
                          const vox::Vector3D& oldVal) {
        if (sdf < 0.05) {
            return vox::Vector3D(0.5, oldVal.y, oldVal.z);
        } else {
            return vox::Vector3D(oldVal);
        }
    });
    
    auto grids = solver->vdbSystemData();
    vox::Size3 resolution = grids->resolution();
    vox::Array2<double> output(resolution.x, resolution.y);
    auto density = solver->smokeDensity();
    char filename[256];
    
    for (vox::Frame frame(0, 1.0 / 60.0); frame.index < 240; ++frame) {
        solver->update(frame);
        
        output.set(0.0);
        density->forEachDataPointIndex(
                                       [&] (const openvdb::Coord& coord) {
            output(coord.x(), coord.y()) += (*density).getGrid()->tree().getValue(coord);
        });
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#grid2,%04d.npy",
                 frame.index);
        saveData(output.constAccessor(), filename);
    }
}
JET_END_TEST_F

JET_BEGIN_TEST_F(GridSmokeSolver3, RisingWithCollider) {
    //
    // This is a replica of smoke_sim example 1.
    //
    
    uint resolutionX = 50;
    
    // Build solver
    auto solver = GridSmokeSolver3::builder()
    .withResolution({resolutionX, 2 * resolutionX, resolutionX})
    .withDomainSizeX(1.0)
    .makeShared();
    
    auto grids = solver->vdbSystemData();
    vox::BoundingBox3D domain = grids->boundingBox();
    
    // Build emitter
    auto box = vox::Box3::builder()
    .withLowerCorner({0.45, -1, 0.45})
    .withUpperCorner({0.55, 0.05, 0.55})
    .makeShared();
    
    auto particle_emitter = VolumeParticleEmitter3::builder()
    .withSurface(box)
    .withSpacing(1.0 / 1024.0)
    .withIsOneShot(true)
    .makeShared();
    solver->setParticleEmitter(particle_emitter);
    
    auto emitter = GridEmitter3::builder()
    .withSourceRegion(box)
    .withIsOneShot(false)
    .makeShared();
    
    solver->setEmitter(emitter);
    emitter->addStepFunctionTarget(solver->smokeDensity(), 0, 1);
    emitter->addStepFunctionTarget(solver->temperature(), 0, 1);
    
    // Build collider
    auto sphere = vox::Sphere3::builder()
    .withCenter({0.5, 0.3, 0.5})
    .withRadius(0.075 * domain.width())
    .makeShared();
    
    auto collider = vox::RigidBodyCollider3::builder()
    .withSurface(sphere)
    .makeShared();
    
    solver->setCollider(collider);
    
    vox::Size3 resolution = grids->resolution();
    vox::Array2<double> output(resolution.x, resolution.y);
    auto density = solver->smokeDensity();
    char filename[256];
    
    for (vox::Frame frame(0, 1.0 / 60.0); frame.index < 240; ++frame) {
        solver->update(frame);
        
        output.set(0.0);
        density->forEachDataPointIndex(
                                       [&] (const openvdb::Coord& coord) {
            output(coord.x(), coord.y()) += (*density).getGrid()->tree().getValue(coord);
        });
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#grid2,%04d.npy",
                 frame.index);
        saveData(output.constAccessor(), filename);
    }
}
JET_END_TEST_F

JET_BEGIN_TEST_F(GridSmokeSolver3, RisingWithColliderLinear) {
    //
    // This is a replica of smoke_sim example 2.
    //
    
    uint resolutionX = 50;
    
    // Build solver
    auto solver = GridSmokeSolver3::builder()
    .withResolution({resolutionX, 2 * resolutionX, resolutionX})
    .withDomainSizeX(1.0)
    .makeShared();
        
    auto grids = solver->vdbSystemData();
    vox::BoundingBox3D domain = grids->boundingBox();
    
    // Build emitter
    auto box = vox::Box3::builder()
    .withLowerCorner({0.45, -1, 0.45})
    .withUpperCorner({0.55, 0.05, 0.55})
    .makeShared();
    
    auto particle_emitter = VolumeParticleEmitter3::builder()
    .withSurface(box)
    .withSpacing(1.0 / 1024.0)
    .withIsOneShot(true)
    .makeShared();
    solver->setParticleEmitter(particle_emitter);
    
    auto emitter = GridEmitter3::builder()
    .withSourceRegion(box)
    .withIsOneShot(false)
    .makeShared();
    
    solver->setEmitter(emitter);
    emitter->addStepFunctionTarget(solver->smokeDensity(), 0, 1);
    emitter->addStepFunctionTarget(solver->temperature(), 0, 1);
    
    // Build collider
    auto sphere = vox::Sphere3::builder()
    .withCenter({0.5, 0.3, 0.5})
    .withRadius(0.075 * domain.width())
    .makeShared();
    
    auto collider = vox::RigidBodyCollider3::builder()
    .withSurface(sphere)
    .makeShared();
    
    solver->setCollider(collider);
    
    vox::Size3 resolution = grids->resolution();
    vox::Array2<double> output(resolution.x, resolution.y);
    auto density = solver->smokeDensity();
    char filename[256];
    
    for (vox::Frame frame(0, 1.0 / 60.0); frame.index < 240; ++frame) {
        solver->update(frame);
        
        output.set(0.0);
        density->forEachDataPointIndex(
                                       [&] (const openvdb::Coord& coord) {
            output(coord.x(), coord.y()) += (*density).getGrid()->tree().getValue(coord);
        });
        snprintf(
                 filename,
                 sizeof(filename),
                 "data.#grid2,%04d.npy",
                 frame.index);
        saveData(output.constAccessor(), filename);
    }
}
JET_END_TEST_F

*/
