//
//  vdb_fmm_level_set_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/9.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_FMM_LEVEL_SET_SOLVER3_H_
#define INCLUDE_VDB_FMM_LEVEL_SET_SOLVER3_H_

#include "../src.common/array_accessor3.h"
#include "vdb_level_set_solver3.hpp"
#include <memory>

namespace vdb {

//!
//! \brief Three-dimensional fast marching method (FMM) implementation.
//!
//! This class implements 3-D FMM. First-order upwind-style differencing is used
//! to solve the PDE.
//!
//! \see https://math.berkeley.edu/~sethian/2006/Explanations/fast_marching_explain.html
//! \see Sethian, James A. "A fast marching level set method for monotonically
//!     advancing fronts." Proceedings of the National Academy of Sciences 93.4
//!     (1996): 1591-1595.
//!
class FmmLevelSetSolver3 final : public LevelSetSolver3 {
public:
    //! Default constructor.
    FmmLevelSetSolver3();
    
    //!
    //! Reinitializes given scalar field to signed-distance field.
    //!
    //! \param inputSdf Input signed-distance field which can be distorted.
    //! \param maxDistance Max range of reinitialization.
    //! \param outputSdf Output signed-distance field.
    //!
    void reinitialize(
                      const ScalarGrid3& inputSdf,
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
    void extrapolate(
                     const ScalarGrid3& input,
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
    void extrapolate(
                     const CollocatedVectorGrid3& input,
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
    void extrapolate(
                     const FaceCenteredGrid3& input,
                     const vox::ScalarField3& sdf,
                     double maxDistance,
                     FaceCenteredGrid3* output) override;
    
private:
    void extrapolate(
                     const openvdb::DoubleGrid::Ptr& input,
                     const vox::ConstArrayAccessor3<double>& sdf,
                     const vox::Vector3D& gridSpacing,
                     double maxDistance,
                     openvdb::DoubleGrid::Ptr output);
};

//! Shared pointer type for the FmmLevelSetSolver3.
typedef std::shared_ptr<FmmLevelSetSolver3> FmmLevelSetSolver3Ptr;

}  // namespace vox

#endif  // INCLUDE_VDB_FMM_LEVEL_SET_SOLVER3_H_

