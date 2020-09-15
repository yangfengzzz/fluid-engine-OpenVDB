//
//  vdb_advection_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/21.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_ADVECTION_SOLVER3_H_
#define INCLUDE_VDB_ADVECTION_SOLVER3_H_

#include "vdb_collocated_vector_grid3.h"
#include "../src.common/constant_scalar_field3.h"
#include "../src.common/constants.h"
#include "vdb_scalar_grid3.h"
#include <limits>
#include <memory>

namespace vdb {

//!
//! \brief Abstract based class for 3-D grid-based advection solver.
//!
//! The implementation of this abstract base class should solve 3-D advection
//! equation for scalar and vector fields.
//!
class AdvectionSolver3 {
public:
    AdvectionSolver3();
    
    virtual ~AdvectionSolver3();
    
    //!
    //! \brief Solves advection equation for given scalar grid.
    //!
    //! The implementation of this virtual function should solve advection
    //! equation for given scalar field \p input and underlying vector field
    //! \p flow that carries the input field. The solution after solving the
    //! equation for given time-step \p dt should be stored in scalar field
    //! \p output. The boundary interface is given by a signed-distance field.
    //! The field is negative inside the boundary. By default, a constant field
    //! with max double value (kMaxD or std::numeric_limists<double>::max())
    //! is used, meaning no boundary.
    //!
    //! \param input Input scalar grid.
    //! \param flow Vector field that advects the input field.
    //! \param dt Time-step for the advection.
    //! \param output Output scalar grid.
    //! \param boundarySdf Boundary interface defined by signed-distance
    //!     field.
    //!
    virtual void advect(
                        const ScalarGrid3& input,
                        const vox::VectorField3& flow,
                        double dt,
                        ScalarGrid3* output,
                        const vox::ScalarField3& boundarySdf
                        = vox::ConstantScalarField3(vox::kMaxD)) = 0;
    
    //!
    //! \brief Solves advection equation for given collocated vector grid.
    //!
    //! The implementation of this virtual function should solve advection
    //! equation for given collocated vector grid \p input and underlying vector
    //! field \p flow that carries the input field. The solution after solving
    //! the equation for given time-step \p dt should be stored in vector field
    //! \p output. The boundary interface is given by a signed-distance field.
    //! The field is negative inside the boundary. By default, a constant field
    //! with max double value (kMaxD or std::numeric_limists<double>::max())
    //! is used, meaning no boundary.
    //!
    //! \param input Input vector grid.
    //! \param flow Vector field that advects the input field.
    //! \param dt Time-step for the advection.
    //! \param output Output vector grid.
    //! \param boundarySdf Boundary interface defined by signed-distance
    //!     field.
    //!
    virtual void advect(
                        const CollocatedVectorGrid3& input,
                        const vox::VectorField3& flow,
                        double dt,
                        CollocatedVectorGrid3* output,
                        const vox::ScalarField3& boundarySdf
                        = vox::ConstantScalarField3(vox::kMaxD));
};

//! Shared pointer type for the 3-D advection solver.
typedef std::shared_ptr<AdvectionSolver3> AdvectionSolver3Ptr;

}  // namespace vox

#endif  // INCLUDE_VDB_ADVECTION_SOLVER3_H_
