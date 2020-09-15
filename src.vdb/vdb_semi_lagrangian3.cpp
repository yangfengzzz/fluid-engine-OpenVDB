//
//  vdb_semi_lagrangian3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/21.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_samplers3.h"
#include "../src.common/parallel.h"
#include "vdb_semi_lagrangian3.hpp"
#include <algorithm>

using namespace vdb;

SemiLagrangian3::SemiLagrangian3() {
}

SemiLagrangian3::~SemiLagrangian3() {
}

void SemiLagrangian3::advect(
                             const ScalarGrid3& input,
                             const vox::VectorField3& flow,
                             double dt,
                             ScalarGrid3* output,
                             const vox::ScalarField3& boundarySdf) {
    auto outputDataPos = output->dataPosition();
    auto inputSamplerFunc = getScalarSamplerFunc(input);
    auto inputDataPos = input.dataPosition();
    
    double h = vox::min3(
                         output->gridSpacing().x,
                         output->gridSpacing().y,
                         output->gridSpacing().z);
    
    output->parallelForEachDataPointIndex([&](const openvdb::Coord& coord) {
        if (boundarySdf.sample(inputDataPos(coord)) > 0.0) {
            vox::Vector3D pt = backTrace(
                                         flow, dt, h,
                                         outputDataPos(coord), boundarySdf);
            output->getGrid()->tree().setValueOn(coord, inputSamplerFunc(pt));
        }
    });
}


void SemiLagrangian3::advect(
                             const CollocatedVectorGrid3& input,
                             const vox::VectorField3& flow,
                             double dt,
                             CollocatedVectorGrid3* output,
                             const vox::ScalarField3& boundarySdf) {
    auto inputSamplerFunc = getVectorSamplerFunc(input);
    
    double h = vox::min3(
                         output->gridSpacing().x,
                         output->gridSpacing().y,
                         output->gridSpacing().z);
    
    auto outputDataPos = output->dataPosition();
    auto inputDataPos = input.dataPosition();
    
    output->parallelForEachDataPointIndex([&](const openvdb::Coord& coord) {
        if (boundarySdf.sample(inputDataPos(coord)) > 0.0) {
            vox::Vector3D pt = backTrace(
                                         flow, dt, h,
                                         outputDataPos(coord), boundarySdf);
            vox::Vector3D val = inputSamplerFunc(pt);
            output->getGrid()->tree().setValueOn(coord,
                                                 openvdb::Vec3d(val.x,val.y,val.z));
        }
    });
}

vox::Vector3D SemiLagrangian3::backTrace(
                                         const vox::VectorField3& flow,
                                         double dt,
                                         double h,
                                         const vox::Vector3D& startPt,
                                         const vox::ScalarField3& boundarySdf) {
    
    double remainingT = dt;
    vox::Vector3D pt0 = startPt;
    vox::Vector3D pt1 = startPt;
    
    while (remainingT > vox::kEpsilonD) {
        // Adaptive time-stepping
        vox::Vector3D vel0 = flow.sample(pt0);
        double numSubSteps
        = std::max(std::ceil(vel0.length() * remainingT / h), 1.0);
        dt = remainingT / numSubSteps;
        
        // Mid-point rule
        vox::Vector3D midPt = pt0 - 0.5 * dt * vel0;
        vox::Vector3D midVel = flow.sample(midPt);
        pt1 = pt0 - dt * midVel;
        
        // Boundary handling
        double phi0 = boundarySdf.sample(pt0);
        double phi1 = boundarySdf.sample(pt1);
        
        if (phi0 * phi1 < 0.0) {
            double w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1));
            pt1 = w * pt0 + (1.0 - w) * pt1;
            break;
        }
        
        remainingT -= dt;
        pt0 = pt1;
    }
    
    return pt1;
}

std::function<double(const vox::Vector3D&)>
SemiLagrangian3::getScalarSamplerFunc(const ScalarGrid3& input) const {
    return input.sampler();
}

std::function<vox::Vector3D(const vox::Vector3D&)>
SemiLagrangian3::getVectorSamplerFunc(
                                      const CollocatedVectorGrid3& input) const {
    return input.sampler();
}
