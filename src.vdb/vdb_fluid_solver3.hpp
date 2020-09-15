//
//  vdb_fluid_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/6.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_FLUID_SOLVER3_H_
#define INCLUDE_VDB_FLUID_SOLVER3_H_

#include "vdb_advection_solver3.hpp"
#include "vdb_cell_centered_scalar_grid3.h"
#include "vdb_face_centered_grid3.h"
#include "../src.common/collider3.h"
#include "vdb_boundary_condition_solver3.hpp"
#include "vdb_diffusion_solver3.hpp"
#include "vdb_pressure_solver3.hpp"
#include "vdb_volume_grid_emitter3.hpp"
#include "vdb_grid_system_data3.h"
#include "../src.common/physics_animation.h"

namespace vdb {

//!
//! \brief Abstract base class for grid-based 3-D fluid solver.
//!
//! This is an abstract base class for grid-based 3-D fluid solver based on
//! Jos Stam's famous 1999 paper - "Stable Fluids".
//!
class FluidSolver3 : public vox::PhysicsAnimation {
public:
    class Builder;
    
    //! Default constructor.
    FluidSolver3();
    
    //! Constructs solver with initial grid size.
    FluidSolver3(const vox::Size3& resolution,
                 const vox::Vector3D& gridSpacing,
                 const vox::Vector3D& gridOrigin);
    
    //! Default destructor.
    virtual ~FluidSolver3();
    
    //! Returns the gravity vector of the system.
    const vox::Vector3D& gravity() const;
    
    //! Sets the gravity of the system.
    void setGravity(const vox::Vector3D& newGravity);
    
    //! Returns the viscosity coefficient.
    double viscosityCoefficient() const;
    
    //!
    //! \brief Sets the viscosity coefficient.
    //!
    //! This function sets the viscosity coefficient. Non-positive input will be
    //! clamped to zero.
    //!
    //! \param[in] newValue The new viscosity coefficient value.
    //!
    void setViscosityCoefficient(double newValue);
    
    //!
    //! \brief Returns the CFL number from the current velocity field for given
    //!     time interval.
    //!
    //! \param[in] timeIntervalInSeconds The time interval in seconds.
    //!
    double cfl(double timeIntervalInSeconds) const;
    
    //! Returns the max allowed CFL number.
    double maxCfl() const;
    
    //! Sets the max allowed CFL number.
    void setMaxCfl(double newCfl);
    
    //! Returns true if the solver is using compressed linear system.
    bool useCompressedLinearSystem() const;
    
    //! Sets whether the solver should use compressed linear system.
    void setUseCompressedLinearSystem(bool onoff);
    
public:
    //! Returns the advection solver instance.
    const AdvectionSolver3Ptr& advectionSolver() const;
    
    //! Sets the advection solver.
    void setAdvectionSolver(const AdvectionSolver3Ptr& newSolver);
    
    //! Returns the diffusion solver instance.
    const DiffusionSolver3Ptr& diffusionSolver() const;
    
    //! Sets the diffusion solver.
    void setDiffusionSolver(const DiffusionSolver3Ptr& newSolver);
    
    //! Returns the pressure solver instance.
    const PressureSolver3Ptr& pressureSolver() const;
    
    //! Sets the pressure solver.
    void setPressureSolver(const PressureSolver3Ptr& newSolver);
    
    //! Returns the closed domain boundary flag.
    int closedDomainBoundaryFlag() const;
    
    //! Sets the closed domain boundary flag.
    void setClosedDomainBoundaryFlag(int flag);
    
    //!
    //! \brief Returns the grid system data.
    //!
    //! This function returns the grid system data. The grid system data stores
    //! the core fluid flow fields such as velocity. By default, the data
    //! instance has velocity field only.
    //!
    //! \see GridSystemData3
    //!
    const GridSystemData3Ptr& vdbSystemData() const;
    
    //!
    //! \brief Resizes grid system data.
    //!
    //! This function resizes grid system data. You can also resize the grid by
    //! calling resize function directly from
    //! GridFluidSolver3::gridSystemData(), but this function provides a
    //! shortcut for the same operation.
    //!
    //! \param[in] newSize        The new size.
    //! \param[in] newGridSpacing The new grid spacing.
    //! \param[in] newGridOrigin  The new grid origin.
    //!
    void resizeGrid(const vox::Size3& newSize,
                    const vox::Vector3D& newGridSpacing,
                    const vox::Vector3D& newGridOrigin);
    
    //!
    //! \brief Returns the resolution of the grid system data.
    //!
    //! This function returns the resolution of the grid system data. This is
    //! equivalent to calling gridSystemData()->resolution(), but provides a
    //! shortcut.
    //!
    vox::Size3 resolution() const;
    
    //!
    //! \brief Returns the grid spacing of the grid system data.
    //!
    //! This function returns the resolution of the grid system data. This is
    //! equivalent to calling gridSystemData()->gridSpacing(), but provides a
    //! shortcut.
    //!
    vox::Vector3D gridSpacing() const;
    
    //!
    //! \brief Returns the origin of the grid system data.
    //!
    //! This function returns the resolution of the grid system data. This is
    //! equivalent to calling gridSystemData()->origin(), but provides a
    //! shortcut.
    //!
    vox::Vector3D gridOrigin() const;
    
    //!
    //! \brief Returns the velocity field.
    //!
    //! This function returns the velocity field from the grid system data.
    //! It is just a shortcut to the most commonly accessed data chunk.
    //!
    const FaceCenteredGrid3Ptr& velocity() const;
    
    //! Returns the collider.
    const vox::Collider3Ptr& collider() const;
    
    //! Sets the collider.
    void setCollider(const vox::Collider3Ptr& newCollider);
    
    //! Returns the emitter.
    const vox::GridEmitter3Ptr& emitter() const;
    
    //! Sets the emitter.
    void setEmitter(const vox::GridEmitter3Ptr& newEmitter);
    
    //! Returns builder fox GridFluidSolver3.
    static Builder builder();
    
protected:
    //! Called when it needs to setup initial condition.
    void onInitialize() override;
    
    //! Called when advancing a single time-step.
    void onAdvanceTimeStep(double timeIntervalInSeconds) override;
    
    //!
    //! \brief Returns the required sub-time-steps for given time interval.
    //!
    //! This function returns the required sub-time-steps for given time
    //! interval based on the max allowed CFL number. If the time interval is
    //! too large so that it makes the CFL number greater than the max value,
    //! This function will return a numebr that is greater than 1.
    //!
    //! \see GridFluidSolver3::maxCfl
    //!
    unsigned int numberOfSubTimeSteps(
                                      double timeIntervalInSeconds) const override;
    
    //! Called at the beginning of a time-step.
    virtual void onBeginAdvanceTimeStep(double timeIntervalInSeconds);
    
    //! Called at the end of a time-step.
    virtual void onEndAdvanceTimeStep(double timeIntervalInSeconds);
    
    //!
    //! \brief Computes the external force terms.
    //!
    //! This function computes the external force applied for given time
    //! interval. By default, it only computes the gravity.
    //!
    //! \see GridFluidSolver3::computeGravity
    //!
    virtual void computeExternalForces(double timeIntervalInSeconds);
    
    //! Computes the viscosity term using the diffusion solver.
    virtual void computeViscosity(double timeIntervalInSeconds);
    
    //! Computes the pressure term using the pressure solver.
    virtual void computePressure(double timeIntervalInSeconds);
    
    //! Computes the advection term using the advection solver.
    virtual void computeAdvection(double timeIntervalInSeconds);
    
    //! Computes the advection term using the advection solver.
    void advect(double timeIntervalInSeconds);
    
    //! Computes the gravity term.
    void computeGravity(double timeIntervalInSeconds);
    
    //!
    //! \brief Applies the boundary condition to the velocity field.
    //!
    //! This function applies the boundary condition to the velocity field by
    //! constraining the flow based on the boundary condition solver.
    //!
    void applyBoundaryCondition();
    
    //! Extrapolates given field into the collider-occupied region.
    void extrapolateIntoCollider(ScalarGrid3* grid);
    
    //! Extrapolates given field into the collider-occupied region.
    void extrapolateIntoCollider(CollocatedVectorGrid3* grid);
    
    //! Extrapolates given field into the collider-occupied region.
    void extrapolateIntoCollider(FaceCenteredGrid3* grid);
    
protected:
    //!
    //! \breif Returns the signed-distance representation of the fluid.
    //!
    //! This function returns the signed-distance representation of the fluid.
    //! Positive sign area is considered to be atmosphere and won't be included
    //! for computing the dynamics. By default, this will return constant scalar
    //! field of -kMaxD, meaning that the entire volume is occupied with fluid.
    //!
    virtual vox::ScalarField3Ptr fluidSdf() const;
    
    //! Returns the signed-distance field representation of the collider.
    vox::ScalarField3Ptr colliderSdf() const;
    
    //! Returns the velocity field of the collider.
    vox::VectorField3Ptr colliderVelocityField() const;
    
protected:
    vox::Vector3D _gravity = vox::Vector3D(0.0, -9.8, 0.0);
    double _viscosityCoefficient = 0.0;
    double _maxCfl = 5.0;
    bool _useCompressedLinearSys = false;
    int _closedDomainBoundaryFlag = vox::kDirectionAll;
    
    GridSystemData3Ptr _grids;
    vox::Collider3Ptr _collider;
    vox::GridEmitter3Ptr _emitter;
    
    AdvectionSolver3Ptr _advectionSolver;
    DiffusionSolver3Ptr _diffusionSolver;
    PressureSolver3Ptr _pressureSolver;
    BoundaryConditionSolver3Ptr _boundaryConditionSolver;
    
    void beginAdvanceTimeStep(double timeIntervalInSeconds);
    
    void endAdvanceTimeStep(double timeIntervalInSeconds);
    
    void updateCollider(double timeIntervalInSeconds);
    
    void updateEmitter(double timeIntervalInSeconds);
};

//! Shared pointer type for the GridFluidSolver3.
typedef std::shared_ptr<FluidSolver3> FluidSolver3Ptr;


//!
//! \brief Base class for grid-based fluid solver builder.
//!
template <typename DerivedBuilder>
class GridFluidSolverBuilderBase3 {
public:
    //! Returns builder with grid resolution.
    DerivedBuilder& withResolution(const vox::Size3& resolution);
    
    //! Returns builder with grid spacing.
    DerivedBuilder& withGridSpacing(const vox::Vector3D& gridSpacing);
    
    //! Returns builder with grid spacing.
    DerivedBuilder& withGridSpacing(double gridSpacing);
    
    //!
    //! \brief Returns builder with domain size in x-direction.
    //!
    //! To build a solver, one can use either grid spacing directly or domain
    //! size in x-direction to set the final grid spacing.
    //!
    DerivedBuilder& withDomainSizeX(double domainSizeX);
    
    //! Returns builder with grid origin
    DerivedBuilder& withOrigin(const vox::Vector3D& gridOrigin);
    
protected:
    vox::Size3 _resolution{1, 1, 1};
    vox::Vector3D _gridSpacing{1, 1, 1};
    vox::Vector3D _gridOrigin{0, 0, 0};
    double _domainSizeX = 1.0;
    bool _useDomainSize = false;
    
    vox::Vector3D getGridSpacing() const;
};

template <typename T>
T& GridFluidSolverBuilderBase3<T>::withResolution(const vox::Size3& resolution) {
    _resolution = resolution;
    return static_cast<T&>(*this);
}

template <typename T>
T& GridFluidSolverBuilderBase3<T>::withGridSpacing(
                                                   const vox::Vector3D& gridSpacing) {
    _gridSpacing = gridSpacing;
    _useDomainSize = false;
    return static_cast<T&>(*this);
}

template <typename T>
T& GridFluidSolverBuilderBase3<T>::withGridSpacing(double gridSpacing) {
    _gridSpacing.x() = gridSpacing;
    _gridSpacing.y() = gridSpacing;
    _gridSpacing.z() = gridSpacing;
    _useDomainSize = false;
    return static_cast<T&>(*this);
}

template <typename T>
T& GridFluidSolverBuilderBase3<T>::withDomainSizeX(double domainSizeX) {
    _domainSizeX = domainSizeX;
    _useDomainSize = true;
    return static_cast<T&>(*this);
}

template <typename T>
T& GridFluidSolverBuilderBase3<T>::withOrigin(const vox::Vector3D& gridOrigin) {
    _gridOrigin = gridOrigin;
    return static_cast<T&>(*this);
}

template <typename T>
vox::Vector3D GridFluidSolverBuilderBase3<T>::getGridSpacing() const {
    vox::Vector3D gridSpacing = _gridSpacing;
    if (_useDomainSize) {
        gridSpacing.set(_domainSizeX / static_cast<double>(_resolution.x));
    }
    return gridSpacing;
}

//!
//! \brief Front-end to create GridFluidSolver3 objects step by step.
//!
class FluidSolver3::Builder final
: public GridFluidSolverBuilderBase3<FluidSolver3::Builder> {
public:
    //! Builds GridFluidSolver3.
    FluidSolver3 build() const;
    
    //! Builds shared pointer of GridFluidSolver3 instance.
    FluidSolver3Ptr makeShared() const {
        return std::make_shared<FluidSolver3>(_resolution,
                                              getGridSpacing(),
                                              _gridOrigin);
    }
};

}  // namespace vox
#endif /* vdb_fluid_solver3_hpp */
