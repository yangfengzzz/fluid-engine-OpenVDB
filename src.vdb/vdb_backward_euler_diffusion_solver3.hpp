//
//  vdb_backward_euler_diffusion_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_BACKWARD_EULER_DIFFUSION_SOLVER3_H_
#define INCLUDE_VDB_BACKWARD_EULER_DIFFUSION_SOLVER3_H_

#include "../src.common/constant_scalar_field3.h"
#include "../src.common/fdm_linear_system_solver3.h"
#include "vdb_diffusion_solver3.hpp"
#include <limits>
#include <memory>

namespace vdb {

//!
//! \brief 3-D grid-based backward Euler diffusion solver.
//!
//! This class implements 3-D grid-based forward Euler diffusion solver using
//! second-order central differencing spatially. Since the method is following
//! the implicit time-integration (i.e. backward Euler), larger time interval or
//! diffusion coefficient can be used without breaking the result. Note, higher
//! values for those parameters will still impact the accuracy of the result.
//! To solve the backward Euler method, a linear system solver is used and
//! incomplete Cholesky conjugate gradient method is used by default.
//!
class BackwardEulerDiffusionSolver3 final : public DiffusionSolver3 {
public:
    enum BoundaryType {
        Dirichlet,
        Neumann
    };
    
    //! Constructs the solver with given boundary type.
    explicit BackwardEulerDiffusionSolver3(
                                              BoundaryType boundaryType
                                              = Neumann);
    
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
    BoundaryType _boundaryType;
    vox::FdmLinearSystem3 _system;
    vox::FdmLinearSystemSolver3Ptr _systemSolver;
    vox::Array3<char> _markers;
    
    void buildMarkers(
                      const vox::Size3& size,
                      const std::function<vox::Vector3D(const openvdb::Coord& coord)>& pos,
                      const vox::ScalarField3& boundarySdf,
                      const vox::ScalarField3& fluidSdf);
    
    void buildMatrix(
                     const vox::Size3& size,
                     const vox::Vector3D& c);
    
    void buildVectors(
                      const openvdb::DoubleGrid::Ptr f,
                      const vox::Size3 size,
                      const vox::Vector3D& c);
    
    void buildVectors(
                      const openvdb::Vec3dGrid::Ptr f,
                      const vox::Size3 size,
                      const vox::Vector3D& c,
                      uint component);
};

//! Shared pointer type for the GridBackwardEulerDiffusionSolver3.
typedef std::shared_ptr<BackwardEulerDiffusionSolver3>
BackwardEulerDiffusionSolver3Ptr;

}  // namespace vox
#endif /* vdb_backward_euler_diffusion_solver3_hpp */
