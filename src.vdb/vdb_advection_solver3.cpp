//
//  vdb_advection_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/21.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_advection_solver3.hpp"

#include <limits>

using namespace vdb;

AdvectionSolver3::AdvectionSolver3() {
}

AdvectionSolver3::~AdvectionSolver3() {
}

void AdvectionSolver3::advect(
                              const CollocatedVectorGrid3& source,
                              const vox::VectorField3& flow,
                              double dt,
                              CollocatedVectorGrid3* target,
                              const vox::ScalarField3& boundarySdf) {
    UNUSED_VARIABLE(source);
    UNUSED_VARIABLE(flow);
    UNUSED_VARIABLE(dt);
    UNUSED_VARIABLE(target);
    UNUSED_VARIABLE(boundarySdf);
}
