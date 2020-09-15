//
//  vdb_boundary_condition_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_BOUNDARY_CONDITION_SOLVER3_H_
#define INCLUDE_VDB_BOUNDARY_CONDITION_SOLVER3_H_

#include "../src.common/collider3.h"
#include "../src.common/constants.h"
#include "vdb_face_centered_grid3.h"
#include "../src.common/scalar_field3.h"

#include <memory>

namespace vdb {

//!
//! \brief Abstract base class for 3-D boundary condition solver for grids.
//!
//! This is a helper class to constrain the 3-D velocity field with given
//! collider object. It also determines whether to open any domain boundaries.
//! To control the friction level, tune the collider parameter.
//!
class BoundaryConditionSolver3 {
public:
    //! Default constructor.
    BoundaryConditionSolver3();
    
    //! Default destructor.
    virtual ~BoundaryConditionSolver3();
    
    //! Returns associated collider.
    const vox::Collider3Ptr& collider() const;
    
    //!
    //! \brief Applies new collider and build the internals.
    //!
    //! This function is called to apply new collider and build the internal
    //! cache. To provide a hint to the cache, info for the expected velocity
    //! grid that will be constrained is provided.
    //!
    //! \param newCollider New collider to apply.
    //!
    void updateCollider(
                        const vox::Collider3Ptr& newCollider,
                        const vox::Size3& gridSize,
                        const vox::Vector3D& gridSpacing,
                        const vox::Vector3D& gridOrigin);
    
    //! Returns the closed domain boundary flag.
    int closedDomainBoundaryFlag() const;
    
    //! Sets the closed domain boundary flag.
    void setClosedDomainBoundaryFlag(int flag);
    
    //!
    //! Constrains the velocity field to conform the collider boundary.
    //!
    //! \param velocity Input and output velocity grid.
    //! \param extrapolationDepth Number of inner-collider grid cells that
    //!     velocity will get extrapolated.
    //!
    virtual void constrainVelocity(
                                   FaceCenteredGrid3* velocity,
                                   unsigned int extrapolationDepth = 5) = 0;
    
    //! Returns the signed distance field of the collider.
    virtual vox::ScalarField3Ptr colliderSdf() const = 0;
    
    //! Returns the velocity field of the collider.
    virtual vox::VectorField3Ptr colliderVelocityField() const = 0;
    
protected:
    //! Invoked when a new collider is set.
    virtual void onColliderUpdated(const vox::Size3& gridSize,
                                   const vox::Vector3D& gridSpacing,
                                   const vox::Vector3D& gridOrigin) = 0;
    
    //! Returns the size of the velocity grid to be constrained.
    const vox::Size3& gridSize() const;
    
    //! Returns the spacing of the velocity grid to be constrained.
    const vox::Vector3D& gridSpacing() const;
    
    //! Returns the origin of the velocity grid to be constrained.
    const vox::Vector3D& gridOrigin() const;
    
private:
    vox::Collider3Ptr _collider;
    vox::Size3 _gridSize;
    vox::Vector3D _gridSpacing;
    vox::Vector3D _gridOrigin;
    int _closedDomainBoundaryFlag = vox::kDirectionAll;
};

//! Shared pointer type for the GridBoundaryConditionSolver3.
typedef std::shared_ptr<BoundaryConditionSolver3>
BoundaryConditionSolver3Ptr;

}  // namespace vox
#endif /* vdb_boundary_condition_solver3_hpp */
