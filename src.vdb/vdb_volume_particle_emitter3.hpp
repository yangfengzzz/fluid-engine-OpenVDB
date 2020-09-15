// Copyright (c) 2019 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_VOLUME_PARTICLE_EMITTER3_H_
#define INCLUDE_VDB_VOLUME_PARTICLE_EMITTER3_H_

#include "../src.common/bounding_box3.h"
#include "../src.common/implicit_surface3.h"
#include "vdb_particle_emitter3.hpp"
#include "../src.common/point_generator3.h"

#include <limits>
#include <memory>
#include <random>

namespace vdb {

//!
//! \brief 3-D volumetric particle emitter.
//!
//! This class emits particles from volumetric geometry.
//!
class VolumeParticleEmitter3 final : public ParticleEmitter3 {
public:
    class Builder;
    
    //!
    //! Constructs an emitter that spawns particles from given implicit surface
    //! which defines the volumetric geometry. Provided bounding box limits
    //! the particle generation region.
    //!
    //! \param[in]  implicitSurface         The implicit surface.
    //! \param[in]  maxRegion               The max region.
    //! \param[in]  spacing                 The spacing between particles.
    //! \param[in]  initialVel              The initial velocity.
    //! \param[in]  linearVel               The linear velocity of the emitter.
    //! \param[in]  angularVel              The angular velocity of the emitter.
    //! \param[in]  maxNumberOfParticles    The max number of particles to be
    //!                                     emitted.
    //! \param[in]  jitter                  The jitter amount between 0 and 1.
    //! \param[in]  isOneShot               True if emitter gets disabled after one shot.
    //! \param[in]  allowOverlapping        True if particles can be overlapped.
    //! \param[in]  seed                    The random seed.
    //!
    VolumeParticleEmitter3(
                           const vox::ImplicitSurface3Ptr& implicitSurface,
                           const vox::BoundingBox3D& maxRegion,
                           double spacing,
                           const vox::Vector3D& initialVel = vox::Vector3D(),
                           const vox::Vector3D& linearVel = vox::Vector3D(),
                           const vox::Vector3D& angularVel = vox::Vector3D(),
                           size_t maxNumberOfParticles = vox::kMaxSize,
                           double jitter = 0.0,
                           bool isOneShot = true,
                           bool allowOverlapping = false,
                           uint32_t seed = 0);
    
    //!
    //! \brief      Sets the point generator.
    //!
    //! This function sets the point generator that defines the pattern of the
    //! point distribution within the volume.
    //!
    //! \param[in]  newPointsGen The new points generator.
    //!
    void setPointGenerator(const vox::PointGenerator3Ptr& newPointsGen);
    
    //! Returns source surface.
    const vox::ImplicitSurface3Ptr& surface() const;
    
    //! Sets the source surface.
    void setSurface(const vox::ImplicitSurface3Ptr& newSurface);
    
    //! Returns max particle gen region.
    const vox::BoundingBox3D& maxRegion() const;
    
    //! Sets the max particle gen region.
    void setMaxRegion(const vox::BoundingBox3D& newBox);
    
    //! Returns jitter amount.
    double jitter() const;
    
    //! Sets jitter amount between 0 and 1.
    void setJitter(double newJitter);
    
    //! Returns true if particles should be emitted just once.
    bool isOneShot() const;
    
    //!
    //! \brief      Sets the flag to true if particles are emitted just once.
    //!
    //! If true is set, the emitter will generate particles only once even after
    //! multiple emit calls. If false, it will keep generating particles from
    //! the volumetric geometry. Default value is true.
    //!
    //! \param[in]  newValue True if particles should be emitted just once.
    //!
    void setIsOneShot(bool newValue);
    
    //! Returns true if particles can be overlapped.
    bool allowOverlapping() const;
    
    //!
    //! \brief      Sets the flag to true if particles can overlap each other.
    //!
    //! If true is set, the emitter will generate particles even if the new
    //! particles can find existing nearby particles within the particle
    //! spacing.
    //!
    //! \param[in]  newValue True if particles can be overlapped.
    //!
    void setAllowOverlapping(bool newValue);
    
    //! Returns max number of particles to be emitted.
    size_t maxNumberOfParticles() const;
    
    //! Sets the max number of particles to be emitted.
    void setMaxNumberOfParticles(size_t newMaxNumberOfParticles);
    
    //! Returns the spacing between particles.
    double spacing() const;
    
    //! Sets the spacing between particles.
    void setSpacing(double newSpacing);
    
    //! Sets the initial velocity of the particles.
    vox::Vector3D initialVelocity() const;
    
    //! Returns the initial velocity of the particles.
    void setInitialVelocity(const vox::Vector3D& newInitialVel);
    
    //! Returns the linear velocity of the emitter.
    vox::Vector3D linearVelocity() const;
    
    //! Sets the linear velocity of the emitter.
    void setLinearVelocity(const vox::Vector3D& newLinearVel);
    
    //! Returns the angular velocity of the emitter.
    vox::Vector3D angularVelocity() const;
    
    //! Sets the linear velocity of the emitter.
    void setAngularVelocity(const vox::Vector3D& newAngularVel);
    
    //! Returns builder fox VolumeParticleEmitter3.
    static Builder builder();
    
private:
    std::mt19937 _rng;
    
    vox::ImplicitSurface3Ptr _implicitSurface;
    vox::BoundingBox3D _bounds;
    double _spacing;
    vox::Vector3D _initialVel;
    vox::Vector3D _linearVel;
    vox::Vector3D _angularVel;
    vox::PointGenerator3Ptr _pointsGen;
    
    size_t _maxNumberOfParticles = vox::kMaxSize;
    size_t _numberOfEmittedParticles = 0;
    
    double _jitter = 0.0;
    bool _isOneShot = true;
    bool _allowOverlapping = false;
    
    //!
    //! \brief      Emits particles to the particle system data.
    //!
    //! \param[in]  currentTimeInSeconds    Current simulation time.
    //! \param[in]  timeIntervalInSeconds   The time-step interval.
    //!
    void onUpdate(
                  double currentTimeInSeconds,
                  double timeIntervalInSeconds) override;
    
    void emit(
              const ParticleSystemData3Ptr& particles,
              vox::Array1<vox::Vector3D>* newPositions,
              vox::Array1<vox::Vector3D>* newVelocities);
    
    double random();
    
    vox::Vector3D velocityAt(const vox::Vector3D& point) const;
};

//! Shared pointer for the VolumeParticleEmitter3 type.
typedef std::shared_ptr<VolumeParticleEmitter3> VolumeParticleEmitter3Ptr;


//!
//! \brief Front-end to create VolumeParticleEmitter3 objects step by step.
//!
class VolumeParticleEmitter3::Builder final {
public:
    //! Returns builder with implicit surface defining volume shape.
    Builder& withImplicitSurface(const vox::ImplicitSurface3Ptr& implicitSurface);
    
    //! Returns builder with surface defining volume shape.
    Builder& withSurface(const vox::Surface3Ptr& surface);
    
    //! Returns builder with max region.
    Builder& withMaxRegion(const vox::BoundingBox3D& bounds);
    
    //! Returns builder with spacing.
    Builder& withSpacing(double spacing);
    
    //! Returns builder with initial velocity.
    Builder& withInitialVelocity(const vox::Vector3D& initialVel);
    
    //! Returns builder with linear velocity.
    Builder& withLinearVelocity(const vox::Vector3D& linearVel);
    
    //! Returns builder with angular velocity.
    Builder& withAngularVelocity(const vox::Vector3D& angularVel);
    
    //! Returns builder with max number of particles.
    Builder& withMaxNumberOfParticles(size_t maxNumberOfParticles);
    
    //! Returns builder with jitter amount.
    Builder& withJitter(double jitter);
    
    //! Returns builder with one-shot flag.
    Builder& withIsOneShot(bool isOneShot);
    
    //! Returns builder with overlapping flag.
    Builder& withAllowOverlapping(bool allowOverlapping);
    
    //! Returns builder with random seed.
    Builder& withRandomSeed(uint32_t seed);
    
    //! Builds VolumeParticleEmitter3.
    VolumeParticleEmitter3 build() const;
    
    //! Builds shared pointer of VolumeParticleEmitter3 instance.
    VolumeParticleEmitter3Ptr makeShared() const;
    
private:
    vox::ImplicitSurface3Ptr _implicitSurface;
    bool _isBoundSet = false;
    vox::BoundingBox3D _bounds;
    double _spacing = 0.1;
    vox::Vector3D _initialVel;
    vox::Vector3D _linearVel;
    vox::Vector3D _angularVel;
    size_t _maxNumberOfParticles = vox::kMaxSize;
    double _jitter = 0.0;
    bool _isOneShot = true;
    bool _allowOverlapping = false;
    uint32_t _seed = 0;
};

}  // namespace vox

#endif  // INCLUDE_VDB_VOLUME_PARTICLE_EMITTER3_H_
