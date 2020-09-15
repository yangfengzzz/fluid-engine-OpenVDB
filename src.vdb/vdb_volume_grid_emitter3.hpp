//
//  vdb_grid_emitter3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/6.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_GRID_EMITTER3_H
#define INCLUDE_VDB_GRID_EMITTER3_H

#include "../src.common/grid_emitter3.h"
#include "vdb_scalar_grid3.h"
#include "vdb_vector_grid3.h"
#include "../src.common/implicit_surface3.h"

#include <tuple>
#include <vector>

namespace vdb {

//!
//! \brief 3-D -based volumetric emitter.
//!
class GridEmitter3 final : public vox::GridEmitter3 {
public:
    class Builder;
    
    //! Maps to a scalar value for given signed-dist, location, and old value.
    typedef std::function<double(double, const vox::Vector3D&, double)> ScalarMapper;
    
    //! Maps to a vector value for given signed-dist, location, and old value.
    typedef std::function<vox::Vector3D(double, const vox::Vector3D&, const vox::Vector3D&)> VectorMapper;
    
    //!
    //! \brief      Constructs an emitter with a source and is-one-shot flag.
    //!
    //! \param[in]  sourceRegion    Emitting region given by the SDF.
    //! \param[in]  isOneShot       True if emitter gets disabled after one shot.
    //!
    explicit GridEmitter3(
                             const vox::ImplicitSurface3Ptr& sourceRegion,
                             bool isOneShot = true);
    
    //! Destructor.
    virtual ~GridEmitter3();
    
public:
    //! Adds signed-distance target to the scalar grid.
    void addSignedDistanceTarget(const ScalarGrid3Ptr& scalarGridTarget);
    
    //!
    //! \brief      Adds step function target to the scalar grid.
    //!
    //! \param[in]  scalarGridTarget The scalar grid target.
    //! \param[in]  minValue         The minimum value of the step function.
    //! \param[in]  maxValue         The maximum value of the step function.
    //!
    void addStepFunctionTarget(
                               const ScalarGrid3Ptr& scalarGridTarget,
                               double minValue,
                               double maxValue);
    
    //!
    //! \brief      Adds a scalar grid target.
    //!
    //! This function adds a custom target to the emitter. The first parameter
    //! defines which grid should it write to. The second parameter,
    //! \p customMapper, defines how to map signed-distance field from the
    //! volume geometry and location of the point to the final scalar value that
    //! is going to be written to the target grid. The third parameter defines
    //! how to blend the old value from the target grid and the new value from
    //! the mapper function.
    //!
    //! \param[in]  scalarGridTarget The scalar grid target
    //! \param[in]  customMapper     The custom mapper.
    //!
    void addTarget(
                   const ScalarGrid3Ptr& scalarGridTarget,
                   const ScalarMapper& customMapper);
    
    //!
    //! \brief      Adds a vector grid target.
    //!
    //! This function adds a custom target to the emitter. The first parameter
    //! defines which grid should it write to. The second parameter,
    //! \p customMapper, defines how to map sigend-distance field from the
    //! volume geometry and location of the point to the final vector value that
    //! is going to be written to the target grid. The third parameter defines
    //! how to blend the old value from the target grid and the new value from
    //! the mapper function.
    //!
    //! \param[in]  vectorGridTarget The vector grid target
    //! \param[in]  customMapper     The custom mapper.
    //!
    void addTarget(
                   const VectorGrid3Ptr& vectorGridTarget,
                   const VectorMapper& customMapper);
    
    //! Returns implicit surface which defines the source region.
    const vox::ImplicitSurface3Ptr& sourceRegion() const;
    
    //! Returns true if this emits only once.
    bool isOneShot() const;
    
    //! Returns builder fox VolumeGridEmitter3.
    static Builder builder();
    
private:
    typedef std::tuple<ScalarGrid3Ptr, ScalarMapper> ScalarTarget;
    typedef std::tuple<VectorGrid3Ptr, VectorMapper> VectorTarget;
    
    vox::ImplicitSurface3Ptr _sourceRegion;
    bool _isOneShot = true;
    bool _hasEmitted = false;
    std::vector<ScalarTarget> _customScalarTargets;
    std::vector<VectorTarget> _customVectorTargets;
    
    void onUpdate(
                  double currentTimeInSeconds,
                  double timeIntervalInSeconds) override;
    
    void emit();
};

//! Shared pointer type for the VolumeGridEmitter3.
typedef std::shared_ptr<GridEmitter3> GridEmitter3Ptr;


//!
//! \brief Front-end to create VolumeGridEmitter3 objects step by step.
//!
class GridEmitter3::Builder final {
public:
    //! Returns builder with surface defining source region.
    Builder& withSourceRegion(const vox::Surface3Ptr& sourceRegion);
    
    //! Returns builder with one-shot flag.
    Builder& withIsOneShot(bool isOneShot);
    
    //! Builds VolumeGridEmitter3.
    GridEmitter3 build() const;
    
    //! Builds shared pointer of VolumeGridEmitter3 instance.
    GridEmitter3Ptr makeShared() const;
    
private:
    vox::ImplicitSurface3Ptr _sourceRegion;
    bool _isOneShot = true;
};

}  // namespace vox

#endif /* INCLUDE__GRID_EMITTER3_H */
