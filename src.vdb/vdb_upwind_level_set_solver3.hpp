//
//  vdb_upwind_level_set_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/9.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_UPWIND_LEVEL_SET_SOLVER3_H_
#define INCLUDE_VDB_UPWIND_LEVEL_SET_SOLVER3_H_

#include "vdb_iterative_level_set_solver3.hpp"

namespace vdb {

//! Three-dimensional first-order upwind-based iterative level set solver.
class UpwindLevelSetSolver3 final : public IterativeLevelSetSolver3 {
public:
    //! Default constructor.
    UpwindLevelSetSolver3();
    
protected:
    //! Computes the derivatives for given grid point.
    void getDerivatives(vox::ConstArrayAccessor3<double> grid,
                        const vox::Vector3D& gridSpacing,
                        uint i,
                        uint j,
                        uint k,
                        std::array<double, 2>* dx,
                        std::array<double, 2>* dy,
                        std::array<double, 2>* dz) const override;
};

typedef std::shared_ptr<UpwindLevelSetSolver3> UpwindLevelSetSolver3Ptr;

}  // namespace vox

#endif  // INCLUDE_VDB_UPWIND_LEVEL_SET_SOLVER3_H_
