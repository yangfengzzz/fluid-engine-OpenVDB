//
//  vdb_diffusion_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_DIFFUSION_SOLVER3_H_
#define INCLUDE_VDB_DIFFUSION_SOLVER3_H_

#include "vdb_collocated_vector_grid3.h"
#include "../src.common/constant_scalar_field3.h"
#include "../src.common/constants.h"
#include "vdb_face_centered_grid3.h"
#include "vdb_scalar_grid3.h"
#include <limits>
#include <memory>

namespace vdb {

//!
//! \brief Abstract base class for 3-D grid-based diffusion equation solver.
//!
//! This class provides functions to solve the diffusion equation for different
//! types of fields. The target equation can be written as
//! \f$\frac{\partial f}{\partial t} = \mu\nabla^2 f\f$ where \f$\mu\f$ is
//! the diffusion coefficient. The field \f$f\f$ can be either scalar or vector
//! field.
//!

class DiffusionSolver3 {
public:
    //! Default constructor.
    DiffusionSolver3();
    
    //! Default destructor.
    virtual ~DiffusionSolver3();
    
    //!
    //! Solves diffusion equation for a scalar field.
    //!
    //! \param source Input scalar field.
    //! \param diffusionCoefficient Amount of diffusion.
    //! \param timeIntervalInSeconds Small time-interval that diffusion occur.
    //! \param dest Output scalar field.
    //! \param boundarySdf Shape of the solid boundary that is empty by default.
    //! \param fluidSdf Shape of the fluid boundary that is full by default.
    //!
    virtual void solve(
                       const ScalarGrid3& source,
                       double diffusionCoefficient,
                       double timeIntervalInSeconds,
                       ScalarGrid3* dest,
                       const vox::ScalarField3& boundarySdf = vox::ConstantScalarField3(vox::kMaxD),
                       const vox::ScalarField3& fluidSdf = vox::ConstantScalarField3(-vox::kMaxD)) = 0;
    
    //!
    //! Solves diffusion equation for a collocated vector field.
    //!
    //! \param source Input collocated vector field.
    //! \param diffusionCoefficient Amount of diffusion.
    //! \param timeIntervalInSeconds Small time-interval that diffusion occur.
    //! \param dest Output collocated vector field.
    //! \param boundarySdf Shape of the solid boundary that is empty by default.
    //! \param fluidSdf Shape of the fluid boundary that is full by default.
    //!
    virtual void solve(
                       const CollocatedVectorGrid3& source,
                       double diffusionCoefficient,
                       double timeIntervalInSeconds,
                       CollocatedVectorGrid3* dest,
                       const vox::ScalarField3& boundarySdf = vox::ConstantScalarField3(vox::kMaxD),
                       const vox::ScalarField3& fluidSdf = vox::ConstantScalarField3(-vox::kMaxD)) = 0;
    
    //!
    //! Solves diffusion equation for a face-centered vector field.
    //!
    //! \param source Input face-centered vector field.
    //! \param diffusionCoefficient Amount of diffusion.
    //! \param timeIntervalInSeconds Small time-interval that diffusion occur.
    //! \param dest Output face-centered vector field.
    //! \param boundarySdf Shape of the solid boundary that is empty by default.
    //! \param fluidSdf Shape of the fluid boundary that is full by default.
    //!
    virtual void solve(
                       const FaceCenteredGrid3& source,
                       double diffusionCoefficient,
                       double timeIntervalInSeconds,
                       FaceCenteredGrid3* dest,
                       const vox::ScalarField3& boundarySdf = vox::ConstantScalarField3(vox::kMaxD),
                       const vox::ScalarField3& fluidSdf = vox::ConstantScalarField3(-vox::kMaxD)) = 0;
    
};

//! Shared pointer type for the GridDiffusionSolver3.
typedef std::shared_ptr<DiffusionSolver3> DiffusionSolver3Ptr;

}  // namespace

#endif /* INCLUDE__DIFFUSION_SOLVER3_H_ */
