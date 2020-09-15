//
//  vdb_cubic_semi_lagrangian3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/21.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_CUBIC_SEMI_LAGRANGIAN3_H_
#define INCLUDE_VDB_CUBIC_SEMI_LAGRANGIAN3_H_

#include "vdb_semi_lagrangian3.hpp"

namespace vdb {

//!
//! \brief Implementation of 3-D cubic semi-Lagrangian advection solver.
//!
//! This class implements 3rd-order cubic 3-D semi-Lagrangian advection solver.
//!
class CubicSemiLagrangian3 final : public SemiLagrangian3 {
public:
    CubicSemiLagrangian3();
    
protected:
    //!
    //! \brief Returns spatial interpolation function object for given scalar
    //! grid.
    //!
    //! This function overrides the original function with cubic interpolation.
    //!
    std::function<double(const vox::Vector3D&)>
    getScalarSamplerFunc(
                         const ScalarGrid3& source) const override;
    
    //!
    //! \brief Returns spatial interpolation function object for given
    //! collocated vector grid.
    //!
    //! This function overrides the original function with cubic interpolation.
    //!
    std::function<vox::Vector3D(const vox::Vector3D&)>
    getVectorSamplerFunc(
                         const CollocatedVectorGrid3& source) const override;
    
};

typedef std::shared_ptr<CubicSemiLagrangian3> CubicSemiLagrangian3Ptr;

}  // namespace vox

#endif  // INCLUDE_VDB_CUBIC_SEMI_LAGRANGIAN3_H_

