//
//  vdb_fractional_boundary_condition_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_FRACTIONAL_BOUNDARY_CONDITION_SOLVER3_H_
#define INCLUDE_VDB_FRACTIONAL_BOUNDARY_CONDITION_SOLVER3_H_

#include "vdb_cell_centered_scalar_grid3.h"
#include "../src.common/custom_vector_field3.h"
#include "vdb_boundary_condition_solver3.hpp"

#include <memory>

namespace vdb {

//!
//! \brief Fractional 3-D boundary condition solver for grids.
//!
//! This class constrains the velocity field by projecting the flow to the
//! signed-distance field representation of the collider. This implementation
//! should pair up with GridFractionalSinglePhasePressureSolver3 to provide
//! sub-grid resolutional velocity projection.
//!
class FractionalBoundaryConditionSolver3
: public BoundaryConditionSolver3 {
public:
    //! Default constructor.
    FractionalBoundaryConditionSolver3();
    
    //! Default destructor.
    virtual ~FractionalBoundaryConditionSolver3();
    
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
    
    //! Returns the signed distance field of the collider.
    vox::ScalarField3Ptr colliderSdf() const override;
    
    //! Returns the velocity field of the collider.
    vox::VectorField3Ptr colliderVelocityField() const override;
    
protected:
    //! Invoked when a new collider is set.
    void onColliderUpdated(const vox::Size3& gridSize,
                           const vox::Vector3D& gridSpacing,
                           const vox::Vector3D& gridOrigin) override;
    
private:
    CellCenteredScalarGrid3Ptr _colliderSdf;
    vox::CustomVectorField3Ptr _colliderVel;
};

//! Shared pointer type for the GridFractionalBoundaryConditionSolver3.
typedef std::shared_ptr<FractionalBoundaryConditionSolver3>
FractionalBoundaryConditionSolver3Ptr;

}  // namespace vox
#endif /* vdb_fractional_boundary_condition_solver3_hpp */
