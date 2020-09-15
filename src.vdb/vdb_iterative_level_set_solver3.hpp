//
//  vdb_iterative_level_set_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/9.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_ITERATIVE_LEVEL_SET_SOLVER3_H_
#define INCLUDE_VDB_ITERATIVE_LEVEL_SET_SOLVER3_H_

#include "../src.common/array_accessor3.h"
#include "vdb_level_set_solver3.hpp"

namespace vdb {

//!
//! \brief Abstract base class for 3-D PDE-based iterative level set solver.
//!
//! This class provides infrastructure for 3-D PDE-based iterative level set
//! solver. Internally, the class implements upwind-style wave propagation and
//! the inheriting classes must provide a way to compute the derivatives for
//! given grid points.
//!
//! \see Osher, Stanley, and Ronald Fedkiw. Level set methods and dynamic
//!     implicit surfaces. Vol. 153. Springer Science & Business Media, 2006.
//!
class IterativeLevelSetSolver3 : public LevelSetSolver3 {
public:
    //! Default constructor.
    IterativeLevelSetSolver3();
    
    //! Default destructor.
    virtual ~IterativeLevelSetSolver3();
    
    //!
    //! Reinitializes given scalar field to signed-distance field.
    //!
    //! \param inputSdf Input signed-distance field which can be distorted.
    //! \param maxDistance Max range of reinitialization.
    //! \param outputSdf Output signed-distance field.
    //!
    void reinitialize(const ScalarGrid3& inputSdf,
                      double maxDistance,
                      ScalarGrid3* outputSdf) override;
    
    //!
    //! Extrapolates given scalar field from negative to positive SDF region.
    //!
    //! \param input Input scalar field to be extrapolated.
    //! \param sdf Reference signed-distance field.
    //! \param maxDistance Max range of extrapolation.
    //! \param output Output scalar field.
    //!
    void extrapolate(const ScalarGrid3& input,
                     const vox::ScalarField3& sdf,
                     double maxDistance,
                     ScalarGrid3* output) override;
    
    //!
    //! Extrapolates given collocated vector field from negative to positive SDF
    //! region.
    //!
    //! \param input Input collocated vector field to be extrapolated.
    //! \param sdf Reference signed-distance field.
    //! \param maxDistance Max range of extrapolation.
    //! \param output Output collocated vector field.
    //!
    void extrapolate(const CollocatedVectorGrid3& input,
                     const vox::ScalarField3& sdf,
                     double maxDistance,
                     CollocatedVectorGrid3* output) override;
    
    //!
    //! Extrapolates given face-centered vector field from negative to positive
    //! SDF region.
    //!
    //! \param input Input face-centered field to be extrapolated.
    //! \param sdf Reference signed-distance field.
    //! \param maxDistance Max range of extrapolation.
    //! \param output Output face-centered vector field.
    //!
    void extrapolate(const FaceCenteredGrid3& input,
                     const vox::ScalarField3& sdf,
                     double maxDistance,
                     FaceCenteredGrid3* output) override;
    
    //! Returns the maximum CFL limit.
    double maxCfl() const;
    
    //!
    //! \brief Sets the maximum CFL limit.
    //!
    //! This function sets the maximum CFL limit for the internal upwind-style
    //! PDE calculation. The negative input will be clamped to 0.
    //!
    void setMaxCfl(double newMaxCfl);
    
protected:
    //! Computes the derivatives for given grid point.
    virtual void getDerivatives(vox::ConstArrayAccessor3<double> grid,
                                const vox::Vector3D& gridSpacing,
                                uint i,
                                uint j,
                                uint k,
                                std::array<double, 2>* dx,
                                std::array<double, 2>* dy,
                                std::array<double, 2>* dz) const = 0;
    
private:
    double _maxCfl = 0.5;
    
    void extrapolate(const openvdb::DoubleGrid::Ptr& input,
                     const vox::ConstArrayAccessor3<double>& sdf,
                     const vox::Vector3D& gridSpacing,
                     double maxDistance,
                     openvdb::DoubleGrid::Ptr output);
    
    static unsigned int distanceToNumberOfIterations(double distance,
                                                     double dtau);
    
    static double sign(const vox::ConstArrayAccessor3<double>& sdf,
                       const vox::Vector3D& gridSpacing,
                       uint i,
                       uint j,
                       uint k);
    
    double pseudoTimeStep(const vox::ConstArrayAccessor3<double>& sdf,
                          const vox::Vector3D& gridSpacing);
};

typedef std::shared_ptr<IterativeLevelSetSolver3> IterativeLevelSetSolver3Ptr;

}  // namespace vox

#endif  // INCLUDE_VDB_ITERATIVE_LEVEL_SET_SOLVER3_H_

