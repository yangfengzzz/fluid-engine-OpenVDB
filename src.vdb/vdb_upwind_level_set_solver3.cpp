//
//  vdb_upwind_level_set_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/9.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "../src.common/pde.h"
#include "vdb_upwind_level_set_solver3.hpp"

#include <algorithm>

using namespace vdb;

UpwindLevelSetSolver3::UpwindLevelSetSolver3() {
    setMaxCfl(0.5);
}

void UpwindLevelSetSolver3::getDerivatives(vox::ConstArrayAccessor3<double> grid,
                                           const vox::Vector3D& gridSpacing,
                                           uint i,
                                           uint j,
                                           uint k,
                                           std::array<double, 2>* dx,
                                           std::array<double, 2>* dy,
                                           std::array<double, 2>* dz) const {
    double D0[3];
    vox::Size3 size = grid.size();
    
    const uint im1 = (i < 1) ? 0 : i - 1;
    const uint ip1 = std::min(i + 1, (uint)size.x - 1);
    const uint jm1 = (j < 1) ? 0 : j - 1;
    const uint jp1 = std::min(j + 1, (uint)size.y - 1);
    const uint km1 = (k < 1) ? 0 : k - 1;
    const uint kp1 = std::min(k + 1, (uint)size.z - 1);
    
    D0[0] = grid(im1, j, k);
    D0[1] = grid(i, j, k);
    D0[2] = grid(ip1, j, k);
    *dx = vox::upwind1(D0, gridSpacing.x);
    
    D0[0] = grid(i, jm1, k);
    D0[1] = grid(i, j, k);
    D0[2] = grid(i, jp1, k);
    *dy = vox::upwind1(D0, gridSpacing.y);
    
    D0[0] = grid(i, j, km1);
    D0[1] = grid(i, j, k);
    D0[2] = grid(i, j, kp1);
    *dz = vox::upwind1(D0, gridSpacing.z);
}
