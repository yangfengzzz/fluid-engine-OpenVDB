//
//  vdb_iterative_level_set_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/9.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "../src.common/array_utils.h"
#include "../src.common/fdm_utils.h"
#include "vdb_iterative_level_set_solver3.hpp"
#include "../src.common/parallel.h"
#include <openvdb/math/Stencils.h>
#include <algorithm>
#include <limits>
#include <utility>  // just make cpplint happy..

using namespace vdb;

IterativeLevelSetSolver3::IterativeLevelSetSolver3() {
}

IterativeLevelSetSolver3::~IterativeLevelSetSolver3() {
}

void IterativeLevelSetSolver3::reinitialize(
                                            const ScalarGrid3& inputSdf,
                                            double maxDistance,
                                            ScalarGrid3* outputSdf) {
    const vox::Size3 size = inputSdf.dataSize();
    const vox::Vector3D gridSpacing = inputSdf.gridSpacing();
    
    JET_THROW_INVALID_ARG_IF(!inputSdf.hasSameShape(*outputSdf));
    
    vox::Array3<double> output(size);
    vox::ArrayAccessor3<double> outputAcc = output.accessor();
    outputAcc.parallelForEachIndex([&](uint i, uint j, uint k){
        outputAcc(i, j, k) = inputSdf.getGrid()->tree().getValue(openvdb::Coord(i, j, k));
    });
    
    const double dtau = pseudoTimeStep(
                                       outputAcc, gridSpacing);
    const unsigned int numberOfIterations
    = distanceToNumberOfIterations(maxDistance, dtau);
    
    vox::Array3<double> temp(size);
    vox::ArrayAccessor3<double> tempAcc = temp.accessor();
    
    LOG(INFO) << "Reinitializing with pseudoTimeStep: " << dtau
    << " numberOfIterations: " << numberOfIterations;
    
    for (unsigned int n = 0; n < numberOfIterations; ++n) {
        output.parallelForEachIndex([&](uint i, uint j, uint k) {
            double s = sign(output, gridSpacing, i, j, k);
            
            std::array<double, 2> dx, dy, dz;
            
            getDerivatives(outputAcc, gridSpacing,
                           i, j, k,
                           &dx, &dy, &dz);
            
            // Explicit Euler step
            double val = outputAcc(i, j, k)
            - dtau * std::max(s, 0.0)
            * (std::sqrt(vox::square(std::max(dx[0], 0.0))
                         + vox::square(std::min(dx[1], 0.0))
                         + vox::square(std::max(dy[0], 0.0))
                         + vox::square(std::min(dy[1], 0.0))
                         + vox::square(std::max(dz[0], 0.0))
                         + vox::square(std::min(dz[1], 0.0))) - 1.0)
            - dtau * std::min(s, 0.0)
            * (std::sqrt(vox::square(std::min(dx[0], 0.0))
                         + vox::square(std::max(dx[1], 0.0))
                         + vox::square(std::min(dy[0], 0.0))
                         + vox::square(std::max(dy[1], 0.0))
                         + vox::square(std::min(dz[0], 0.0))
                         + vox::square(std::max(dz[1], 0.0))) - 1.0);
            tempAcc(i, j, k) = val;
        });
        
        std::swap(temp, output);
    }
    
    outputAcc.forEachIndex([&](uint i, uint j, uint k){
        outputSdf->getGrid()->tree().setValueOnly(openvdb::Coord(i, j, k),
                                                  outputAcc(i, j, k));
    });
}

void IterativeLevelSetSolver3::extrapolate(
                                           const ScalarGrid3& input,
                                           const vox::ScalarField3& sdf,
                                           double maxDistance,
                                           ScalarGrid3* output) {
    JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));
    
    vox::Array3<double> sdfGrid(input.dataSize());
    auto pos = input.dataPosition();
    sdfGrid.parallelForEachIndex([&](uint i, uint j, uint k) {
        sdfGrid(i, j, k) = sdf.sample(pos(openvdb::Coord(i, j, k)));
    });
    
    extrapolate(input.getGrid(), sdfGrid,
                input.gridSpacing(), maxDistance,
                output->getGrid());
}

void IterativeLevelSetSolver3::extrapolate(const CollocatedVectorGrid3& input,
                                           const vox::ScalarField3& sdf,
                                           double maxDistance,
                                           CollocatedVectorGrid3* output) {
    JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));
    
    vox::Array3<double> sdfGrid(input.dataSize());
    auto pos = input.dataPosition();
    sdfGrid.parallelForEachIndex([&](uint i, uint j, uint k) {
        sdfGrid(i, j, k) = sdf.sample(pos(openvdb::Coord(i, j, k)));
    });
    
    const vox::Vector3D gridSpacing = input.gridSpacing();
    
    openvdb::DoubleGrid::Ptr u = openvdb::DoubleGrid::create(0.0);
    u->topologyUnion(*input.getGrid());
    openvdb::DoubleGrid::Ptr u0 = openvdb::DoubleGrid::create(0.0);
    u0->topologyUnion(*input.getGrid());
    openvdb::DoubleGrid::Ptr v = openvdb::DoubleGrid::create(0.0);
    v->topologyUnion(*input.getGrid());
    openvdb::DoubleGrid::Ptr v0 = openvdb::DoubleGrid::create(0.0);
    v0->topologyUnion(*input.getGrid());
    openvdb::DoubleGrid::Ptr w = openvdb::DoubleGrid::create(0.0);
    w->topologyUnion(*input.getGrid());
    openvdb::DoubleGrid::Ptr w0 = openvdb::DoubleGrid::create(0.0);
    w0->topologyUnion(*input.getGrid());
    
    sdfGrid.parallelForEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        u->tree().setValueOnly(coord, input.getGrid()->tree().getValue(coord).x());
        v->tree().setValueOnly(coord, input.getGrid()->tree().getValue(coord).y());
        w->tree().setValueOnly(coord, input.getGrid()->tree().getValue(coord).z());
    });
    
    extrapolate(u, sdfGrid, gridSpacing, maxDistance, u0);
    
    extrapolate(v, sdfGrid, gridSpacing, maxDistance, v0);
    
    extrapolate(w, sdfGrid, gridSpacing, maxDistance, w0);
    
    sdfGrid.parallelForEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        output->getGrid()->tree().setValueOnly(coord,
                                               openvdb::Vec3d(u->tree().getValue(coord),
                                                              v->tree().getValue(coord),
                                                              w->tree().getValue(coord)));
    });
}

void IterativeLevelSetSolver3::extrapolate(const FaceCenteredGrid3& input,
                                           const vox::ScalarField3& sdf,
                                           double maxDistance,
                                           FaceCenteredGrid3* output) {
    JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));
    
    const vox::Vector3D gridSpacing = input.gridSpacing();
    
    vox::Array3<double> sdfAtU(input.uSize());
    auto uPos = input.uPosition();
    sdfAtU.parallelForEachIndex([&](uint i, uint j, uint k) {
        sdfAtU(i, j, k) = sdf.sample(uPos(openvdb::Coord(i, j, k)));
    });
    
    extrapolate(input.getUGrid(), sdfAtU,
                gridSpacing, maxDistance,
                output->getUGrid());
    
    auto vPos = input.vPosition();
    vox::Array3<double> sdfAtV(input.vSize());
    sdfAtV.parallelForEachIndex([&](uint i, uint j, uint k) {
        sdfAtV(i, j, k) = sdf.sample(vPos(openvdb::Coord(i, j, k)));
    });
    
    extrapolate(input.getVGrid(), sdfAtV,
                gridSpacing, maxDistance,
                output->getVGrid());
    
    auto wPos = input.wPosition();
    vox::Array3<double> sdfAtW(input.wSize());
    sdfAtW.parallelForEachIndex([&](uint i, uint j, uint k) {
        sdfAtW(i, j, k) = sdf.sample(wPos(openvdb::Coord(i, j, k)));
    });
    
    extrapolate(input.getWGrid(), sdfAtW,
                gridSpacing, maxDistance,
                output->getWGrid());
}

void IterativeLevelSetSolver3::extrapolate(
                                           const openvdb::DoubleGrid::Ptr& input,
                                           const vox::ConstArrayAccessor3<double>& sdf,
                                           const vox::Vector3D& gridSpacing,
                                           double maxDistance,
                                           openvdb::DoubleGrid::Ptr output) {
    const vox::Size3 size = sdf.size();
    
    vox::Array3<double> outputArray(size);
    vox::ArrayAccessor3<double> outputAcc = outputArray.accessor();
    outputAcc.parallelForEachIndex([&](uint i, uint j, uint k){
        outputAcc(i, j, k) = input->tree().getValue(openvdb::Coord(i, j, k));
    });
    
    const double dtau = pseudoTimeStep(sdf, gridSpacing);
    const unsigned int numberOfIterations
    = distanceToNumberOfIterations(maxDistance, dtau);
    
    vox::Array3<double> temp(size);
    vox::ArrayAccessor3<double> tempAcc = temp.accessor();
    
    for (unsigned int n = 0; n < numberOfIterations; ++n) {
        vox::parallelFor(
                         vox::kZeroSize, size.x,
                         vox::kZeroSize, size.y,
                         vox::kZeroSize, size.z,
                         [&](uint i, uint j, uint k) {
            if (sdf(i, j, k) >= 0) {
                std::array<double, 2> dx, dy, dz;
                vox::Vector3D grad = gradient3(sdf, gridSpacing,
                                               i, j, k);
                
                getDerivatives(
                               outputAcc, gridSpacing, i, j, k, &dx, &dy, &dz);
                
                tempAcc(i, j, k) = outputAcc(i, j, k)
                - dtau * (std::max(grad.x, 0.0) * dx[0]
                          + std::min(grad.x, 0.0) * dx[1]
                          + std::max(grad.y, 0.0) * dy[0]
                          + std::min(grad.y, 0.0) * dy[1]
                          + std::max(grad.z, 0.0) * dz[0]
                          + std::min(grad.z, 0.0) * dz[1]);
            } else {
                tempAcc(i, j, k) = outputAcc(i, j, k);
            }
        });
        
        std::swap(tempAcc, outputAcc);
    }
    
    outputAcc.forEachIndex([&](uint i, uint j, uint k){
        output->tree().setValueOnly(openvdb::Coord(i, j, k),
                                    outputAcc(i, j, k));
    });
}

double IterativeLevelSetSolver3::maxCfl() const {
    return _maxCfl;
}

void IterativeLevelSetSolver3::setMaxCfl(double newMaxCfl) {
    _maxCfl = std::max(newMaxCfl, 0.0);
}

unsigned int IterativeLevelSetSolver3::distanceToNumberOfIterations(
                                                                    double distance,
                                                                    double dtau) {
    return static_cast<unsigned int>(std::ceil(distance / dtau));
}

double IterativeLevelSetSolver3::sign(
                                      const vox::ConstArrayAccessor3<double>& sdf,
                                      const vox::Vector3D& gridSpacing,
                                      uint i,
                                      uint j,
                                      uint k) {
    double d = sdf(i, j, k);
    double e = vox::min3(gridSpacing.x, gridSpacing.y, gridSpacing.z);
    return d / std::sqrt(d * d + e * e);
}

double IterativeLevelSetSolver3::pseudoTimeStep(
                                                const vox::ConstArrayAccessor3<double>& sdf,
                                                const vox::Vector3D& gridSpacing) {
    const vox::Size3 size = sdf.size();
    
    const double h = vox::max3(gridSpacing.x, gridSpacing.y, gridSpacing.z);
    
    double maxS = -std::numeric_limits<double>::max();
    double dtau = _maxCfl * h;
    
    for (uint k = 0; k < size.z; ++k) {
        for (uint j = 0; j < size.y; ++j) {
            for (uint i = 0; i < size.x; ++i) {
                double s = sign(sdf, gridSpacing, i, j, k);
                maxS = std::max(s, maxS);
            }
        }
    }
    
    while (dtau * maxS / h > _maxCfl) {
        dtau *= 0.5;
    }
    
    return dtau;
}
