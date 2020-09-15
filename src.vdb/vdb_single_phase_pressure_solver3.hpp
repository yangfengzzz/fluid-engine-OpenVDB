//
//  vdb_single_phase_pressure_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_SINGLE_PHASE_PRESSURE_SOLVER3_H_
#define INCLUDE_VDB_SINGLE_PHASE_PRESSURE_SOLVER3_H_

#include "../src.common/fdm_linear_system_solver3.h"
#include "../src.common/fdm_mg_linear_system3.h"
#include "../src.common/fdm_mg_solver3.h"
#include "vdb_boundary_condition_solver3.hpp"
#include "vdb_pressure_solver3.hpp"

#include <memory>

class PIC_test;

namespace vdb {

//!
//! \brief 3-D single-phase pressure solver.
//!
//! This class implements 3-D single-phase pressure solver. This solver encodes
//! the boundaries like Lego blocks -- if a grid cell center is inside or
//! outside the boundaries, it is either marked as occupied or not.
//! In addition, this class solves single-phase flow, solving the pressure for
//! selected fluid region only and treat other area as an atmosphere region.
//! Thus, the pressure outside the fluid will be set to a constant value and
//! velocity field won't be altered. This solver also computes the fluid
//! boundary in block-like manner; If a grid cell is inside or outside the
//! fluid, it is marked as either fluid or atmosphere. Thus, this solver in
//! general, does not compute subgrid structure.
//!
class SinglePhasePressureSolver3 : public PressureSolver3 {
public:
    //! Default constructor.
    SinglePhasePressureSolver3();
    
    //! Default destructor.
    virtual ~SinglePhasePressureSolver3();
    
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
    //! GridBlockedBoundaryConditionSolver2 will be returned since this pressure
    //! solver encodes boundaries like pixelated Lego blocks.
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
    
    std::vector<vox::Array3<char>> _markers;
    
    void buildMarkers(const vox::Size3& size,
                      const std::function<vox::Vector3D(const openvdb::Coord& coord)>& pos,
                      const vox::ScalarField3& boundarySdf,
                      const vox::ScalarField3& fluidSdf);
    
    void decompressSolution();
    
    virtual void buildSystem(const FaceCenteredGrid3& input,
                             bool useCompressed);
    
    virtual void applyPressureGradient(const FaceCenteredGrid3& input,
                                       FaceCenteredGrid3* output);
    
    friend class ::PIC_test;
};

//! Shared pointer type for the GridSinglePhasePressureSolver3.
typedef std::shared_ptr<SinglePhasePressureSolver3>
SinglePhasePressureSolver3Ptr;

}  // namespace vox
#endif /* vdb_single_phase_pressure_solver3_hpp */
