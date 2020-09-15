//
//  vdb_eno_level_set_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/9.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "../src.common/pde.h"
#include "vdb_eno_level_set_solver3.hpp"

#include <algorithm>

using namespace vdb;

EnoLevelSetSolver3::EnoLevelSetSolver3() {
    setMaxCfl(0.25);
}

void EnoLevelSetSolver3::getDerivatives(vox::ConstArrayAccessor3<double> grid,
                                        const vox::Vector3D& gridSpacing,
                                        uint i,
                                        uint j,
                                        uint k,
                                        std::array<double, 2>* dx,
                                        std::array<double, 2>* dy,
                                        std::array<double, 2>* dz) const{
    double D0[7];
    vox::Size3 size = grid.size();
    
    const uint im3 = (i < 3) ? 0 : i - 3;
    const uint im2 = (i < 2) ? 0 : i - 2;
    const uint im1 = (i < 1) ? 0 : i - 1;
    const uint ip1 = std::min(i + 1, (uint)size.x - 1);
    const uint ip2 = std::min(i + 2, (uint)size.x - 1);
    const uint ip3 = std::min(i + 3, (uint)size.x - 1);
    const uint jm3 = (j < 3) ? 0 : j - 3;
    const uint jm2 = (j < 2) ? 0 : j - 2;
    const uint jm1 = (j < 1) ? 0 : j - 1;
    const uint jp1 = std::min(j + 1, (uint)size.y - 1);
    const uint jp2 = std::min(j + 2, (uint)size.y - 1);
    const uint jp3 = std::min(j + 3, (uint)size.y - 1);
    const uint km3 = (k < 3) ? 0 : k - 3;
    const uint km2 = (k < 2) ? 0 : k - 2;
    const uint km1 = (k < 1) ? 0 : k - 1;
    const uint kp1 = std::min(k + 1, (uint)size.z - 1);
    const uint kp2 = std::min(k + 2, (uint)size.z - 1);
    const uint kp3 = std::min(k + 3, (uint)size.z - 1);
    
    // 3rd-order ENO differencing
    D0[0] = grid(im3, j, k);
    D0[1] = grid(im2, j, k);
    D0[2] = grid(im1, j, k);
    D0[3] = grid(i, j, k);
    D0[4] = grid(ip1, j, k);
    D0[5] = grid(ip2, j, k);
    D0[6] = grid(ip3, j, k);
    *dx = vox::eno3(D0, gridSpacing.x);
    
    D0[0] = grid(i, jm3, k);
    D0[1] = grid(i, jm2, k);
    D0[2] = grid(i, jm1, k);
    D0[3] = grid(i, j, k);
    D0[4] = grid(i, jp1, k);
    D0[5] = grid(i, jp2, k);
    D0[6] = grid(i, jp3, k);
    *dy = vox::eno3(D0, gridSpacing.y);
    
    D0[0] = grid(i, j, km3);
    D0[1] = grid(i, j, km2);
    D0[2] = grid(i, j, km1);
    D0[3] = grid(i, j, k);
    D0[4] = grid(i, j, kp1);
    D0[5] = grid(i, j, kp2);
    D0[6] = grid(i, j, kp3);
    *dz = vox::eno3(D0, gridSpacing.z);
}
