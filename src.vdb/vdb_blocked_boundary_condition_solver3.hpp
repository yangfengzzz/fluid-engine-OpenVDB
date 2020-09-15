//
//  vdb_blocked_boundary_condition_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_BLOCKED_BOUNDARY_CONDITION_SOLVER3_H_
#define INCLUDE_VDB_BLOCKED_BOUNDARY_CONDITION_SOLVER3_H_

#include "vdb_fractional_boundary_condition_solver3.hpp"
#include "../src.common/array3.h"
#include <memory>

namespace vdb {

//!
//! \brief Blocked 3-D boundary condition solver for grids.
//!
//! This class constrains the velocity field by projecting the flow to the
//! blocked representation of the collider. A collider is rasterized into voxels
//! and each face of the collider voxels projects the velocity field onto its
//! face. This implementation should pair up with GridSinglePhasePressureSolver3
//! since the pressure solver assumes blocked boundary representation as well.
//!
class BlockedBoundaryConditionSolver3 final
: public FractionalBoundaryConditionSolver3 {
public:
    //! Default constructor.
    BlockedBoundaryConditionSolver3();
    
    //!
    //! Constrains the velocity field to conform the collider boundary.
    //!
    //! \param velocity Input and output velocity grid.
    //! \param extrapolationDepth Number of inner-collider grid cells that
    //!     velocity will get extrapolated.
    //!
    void constrainVelocity(
                           FaceCenteredGrid3* velocity,
                           unsigned int extrapolationDepth = 5) override;
    
    //! Returns the marker which is 1 if occupied by the collider.
    const vox::Array3<char>& marker() const;
    
protected:
    //! Invoked when a new collider is set.
    void onColliderUpdated(const vox::Size3& gridSize,
                           const vox::Vector3D& gridSpacing,
                           const vox::Vector3D& gridOrigin) override;
    
private:
    vox::Array3<char> _marker;
};

//! Shared pointer type for the GridBlockedBoundaryConditionSolver3.
typedef std::shared_ptr<BlockedBoundaryConditionSolver3>
BlockedBoundaryConditionSolver3Ptr;

}  // namespace vox
#endif /* vdb_blocked_boundary_condition_solver3_hpp */
