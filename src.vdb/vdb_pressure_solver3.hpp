//
//  vdb_pressure_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_PRESSURE_SOLVER3_H_
#define INCLUDE_VDB_PRESSURE_SOLVER3_H_

#include "vdb_collocated_vector_grid3.h"
#include "../src.common/constant_scalar_field3.h"
#include "../src.common/constant_vector_field3.h"
#include "../src.common/constants.h"
#include "vdb_face_centered_grid3.h"
#include "vdb_boundary_condition_solver3.hpp"
#include "vdb_scalar_grid3.h"
#include <memory>

namespace vdb {

//!
//! \brief Abstract base class for 3-D grid-based pressure solver.
//!
//! This class represents a 3-D grid-based pressure solver interface which can
//! be used as a sub-step of GridFluidSolver3. Inheriting classes must implement
//! the core GridPressureSolver3::solve function as well as the helper function
//! GridPressureSolver3::suggestedBoundaryConditionSolver.
//!
class PressureSolver3 {
public:
    //! Default constructor.
    PressureSolver3();
    
    //! Default destructor.
    virtual ~PressureSolver3();
    
    //!
    //! \brief Solves the pressure term and apply it to the velocity field.
    //!
    //! This function takes input velocity field and outputs pressure-applied
    //! velocity field. It also accepts extra arguments such as \p boundarySdf
    //! and \p fluidSdf that represent signed-distance representation of the
    //! boundary and fluid area. The negative region of \p boundarySdf means
    //! it is occupied by solid object. Also, the positive / negative area of
    //! the \p fluidSdf means it is occupied by fluid / atmosphere. If not
    //! specified, constant scalar field with kMaxD will be used for
    //! \p boundarySdf meaning that no boundary at all. Similarly, a constant
    //! field with -kMaxD will be used for \p fluidSdf which means it's fully
    //! occupied with fluid without any atmosphere.
    //!
    //! \param[in]    input                 The input velocity field.
    //! \param[in]    timeIntervalInSeconds The time interval for the sim.
    //! \param[out]   output                The output velocity field.
    //! \param[in]    boundarySdf           The SDF of the boundary.
    //! \param[in]    fluidSdf              The SDF of the fluid/atmosphere.
    //!
    virtual void solve(
                       const FaceCenteredGrid3& input,
                       double timeIntervalInSeconds,
                       FaceCenteredGrid3* output,
                       const vox::ScalarField3& boundarySdf = vox::ConstantScalarField3(vox::kMaxD),
                       const vox::VectorField3& boundaryVelocity = vox::ConstantVectorField3({0, 0, 0}),
                       const vox::ScalarField3& fluidSdf = vox::ConstantScalarField3(-vox::kMaxD),
                       bool useCompressed = false) = 0;
    
    //!
    //! \brief Returns the best boundary condition solver for this solver.
    //!
    //! This function returns the best boundary condition solver that works well
    //! with this pressure solver. Depending on the pressure solver
    //! implementation, different boundary condition solver might be used.
    //!
    virtual BoundaryConditionSolver3Ptr suggestedBoundaryConditionSolver() const = 0;
};

//! Shared pointer type for the GridPressureSolver3.
typedef std::shared_ptr<PressureSolver3> PressureSolver3Ptr;

}  // namespace vox
#endif /* vdb_pressure_solver3_hpp */
