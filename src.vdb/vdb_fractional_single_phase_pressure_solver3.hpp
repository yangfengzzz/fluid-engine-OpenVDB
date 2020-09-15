//
//  vdb_fractional_single_phase_pressure_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef vdb_fractional_single_phase_pressure_solver3_hpp
#define vdb_fractional_single_phase_pressure_solver3_hpp

#include "vdb_cell_centered_scalar_grid3.h"
#include "../src.common/fdm_linear_system_solver3.h"
#include "../src.common/fdm_mg_linear_system3.h"
#include "../src.common/fdm_mg_solver3.h"
#include "vdb_boundary_condition_solver3.hpp"
#include "vdb_pressure_solver3.hpp"
#include "vdb_vertex_centered_scalar_grid3.h"
#include <memory>

class PIC_test;

namespace vdb {

//!
//! \brief 3-D fractional single-phase pressure solver.
//!
//! This class implements 3-D fractional (or variational) single-phase pressure
//! solver. It is called fractional because the solver encodes the boundaries
//! to the grid cells like anti-aliased pixels, meaning that a grid cell will
//! record the partially overlapping boundary as a fractional number.
//! Alternative apporach is to represent boundaries like Lego blocks which is
//! the case for GridSinglePhasePressureSolver2.
//! In addition, this class solves single-phase flow, solving the pressure for
//! selected fluid region only and treat other area as an atmosphere region.
//! Thus, the pressure outside the fluid will be set to a constant value and
//! velocity field won't be altered. This solver also computes the fluid
//! boundary in fractional manner, meaning that the solver tries to capture the
//! subgrid structures. This class uses ghost fluid method for such calculation.
//!
//! \see Batty, Christopher, Florence Bertails, and Robert Bridson.
//!     "A fast variational framework for accurate solid-fluid coupling."
//!     ACM Transactions on Graphics (TOG). Vol. 26. No. 3. ACM, 2007.
//! \see Enright, Doug, et al. "Using the particle level set method and
//!     a second order accurate pressure boundary condition for free surface
//!     flows." ASME/JSME 2003 4th Joint Fluids Summer Engineering Conference.
//!     American Society of Mechanical Engineers, 2003.
//!
class FractionalSinglePhasePressureSolver3 : public PressureSolver3 {
public:
    FractionalSinglePhasePressureSolver3();
    
    virtual ~FractionalSinglePhasePressureSolver3();
    
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
    void solve(
               const FaceCenteredGrid3& input,
               double timeIntervalInSeconds,
               FaceCenteredGrid3* output,
               const vox::ScalarField3& boundarySdf = vox::ConstantScalarField3(vox::kMaxD),
               const vox::VectorField3& boundaryVelocity = vox::ConstantVectorField3({0, 0, 0}),
               const vox::ScalarField3& fluidSdf = vox::ConstantScalarField3(-vox::kMaxD),
               bool useCompressed = false) override;
    
    //!
    //! \brief Returns the best boundary condition solver for this solver.
    //!
    //! This function returns the best boundary condition solver that works well
    //! with this pressure solver. Depending on the pressure solver
    //! implementation, different boundary condition solver might be used. For
    //! this particular class, an instance of
    //! GridFractionalBoundaryConditionSolver3 will be returned.
    //!
    BoundaryConditionSolver3Ptr suggestedBoundaryConditionSolver()
    const override;
    
    //! Returns the linear system solver.
    const vox::FdmLinearSystemSolver3Ptr& linearSystemSolver() const;
    
    //! Sets the linear system solver.
    void setLinearSystemSolver(const vox::FdmLinearSystemSolver3Ptr& solver);
    
    //! Returns the pressure field.
    const vox::FdmVector3& pressure() const;
    
private:
    vox::FdmLinearSystem3 _system;
    vox::FdmCompressedLinearSystem3 _compSystem;
    vox::FdmLinearSystemSolver3Ptr _systemSolver;
    
    vox::FdmMgLinearSystem3 _mgSystem;
    vox::FdmMgSolver3Ptr _mgSystemSolver;
    
    std::vector<vox::Array3<float>> _uWeights;
    std::vector<vox::Array3<float>> _vWeights;
    std::vector<vox::Array3<float>> _wWeights;
    std::vector<vox::Array3<float>> _fluidSdf;
    
    std::function<vox::Vector3D(const vox::Vector3D&)> _boundaryVel;
    
    void buildWeights(const FaceCenteredGrid3& input,
                      const vox::ScalarField3& boundarySdf,
                      const vox::VectorField3& boundaryVelocity,
                      const vox::ScalarField3& fluidSdf);
    
    void decompressSolution();
        
    virtual void buildSystem(const FaceCenteredGrid3& input,
                             bool useCompressed);
    
    virtual void applyPressureGradient(const FaceCenteredGrid3& input,
                                       FaceCenteredGrid3* output);
    
    friend class ::PIC_test;
};

//! Shared pointer type for the GridFractionalSinglePhasePressureSolver3.
typedef std::shared_ptr<FractionalSinglePhasePressureSolver3>
FractionalSinglePhasePressureSolver3Ptr;

}  // namespace vox
#endif /* vdb_fractional_single_phase_pressure_solver3_hpp */
