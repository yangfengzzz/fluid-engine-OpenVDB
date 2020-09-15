//
//  vdb_forward_euler_diffusion_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_FORWARD_EULER_DIFFUSION_SOLVER3_H_
#define INCLUDE_VDB_FORWARD_EULER_DIFFUSION_SOLVER3_H_

#include "../src.common/array3.h"
#include "../src.common/constant_scalar_field3.h"
#include "vdb_diffusion_solver3.hpp"
#include <limits>
#include <memory>

namespace vdb {

//!
//! \brief 3-D grid-based forward Euler diffusion solver.
//!
//! This class implements 3-D grid-based forward Euler diffusion solver using
//! second-order central differencing spatially. Since the method is relying on
//! explicit time-integration (i.e. forward Euler), the diffusion coefficient is
//! limited by the time interval and grid spacing such as:
//! \f$\mu < \frac{h}{12\Delta t} \f$ where \f$\mu\f$, \f$h\f$, and
//! \f$\Delta t\f$ are the diffusion coefficient, grid spacing, and time
//! interval, respectively.
//!
class ForwardEulerDiffusionSolver3 final : public DiffusionSolver3 {
public:
    //! Default constructor.
    ForwardEulerDiffusionSolver3();
    
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
    void solve(
               const ScalarGrid3& source,
               double diffusionCoefficient,
               double timeIntervalInSeconds,
               ScalarGrid3* dest,
               const vox::ScalarField3& boundarySdf
               = vox::ConstantScalarField3(vox::kMaxD),
               const vox::ScalarField3& fluidSdf
               = vox::ConstantScalarField3(-vox::kMaxD)) override;
    
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
    void solve(
               const CollocatedVectorGrid3& source,
               double diffusionCoefficient,
               double timeIntervalInSeconds,
               CollocatedVectorGrid3* dest,
               const vox::ScalarField3& boundarySdf
               = vox::ConstantScalarField3(vox::kMaxD),
               const vox::ScalarField3& fluidSdf
               = vox::ConstantScalarField3(-vox::kMaxD)) override;
    
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
    void solve(
               const FaceCenteredGrid3& source,
               double diffusionCoefficient,
               double timeIntervalInSeconds,
               FaceCenteredGrid3* dest,
               const vox::ScalarField3& boundarySdf
               = vox::ConstantScalarField3(vox::kMaxD),
               const vox::ScalarField3& fluidSdf
               = vox::ConstantScalarField3(-vox::kMaxD)) override;
    
private:
    vox::Array3<char> _markers;
    
    void buildMarkers(
                      const vox::Size3& size,
                      const std::function<vox::Vector3D(const openvdb::Coord& coord)>& pos,
                      const vox::ScalarField3& boundarySdf,
                      const vox::ScalarField3& fluidSdf);
};

//! Shared pointer type for the GridForwardEulerDiffusionSolver3.
typedef std::shared_ptr<ForwardEulerDiffusionSolver3>
ForwardEulerDiffusionSolver3Ptr;

}  // namespace vox
#endif /* INCLUDE_VDB_FORWARD_EULER_DIFFUSION_SOLVER3_H_ */
