//
//  vdb_fmm_level_set_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/9.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "../src.common/fdm_utils.h"
#include "vdb_fmm_level_set_solver3.hpp"
#include "../src.common/level_set_utils.h"
#include <openvdb/math/Stencils.h>
#include <algorithm>
#include <queue>
#include <vector>

using namespace vdb;

static const char kUnknown = 0;
static const char kKnown = 1;
static const char kTrial = 2;

// Find geometric solution near the boundary
inline double solveQuadNearBoundary(const vox::Array3<char>& markers,
                                    openvdb::DoubleGrid::Ptr output,
                                    const vox::Vector3D& gridSpacing,
                                    const vox::Vector3D& invGridSpacingSqr,
                                    double sign,
                                    uint i, uint j, uint k) {
    UNUSED_VARIABLE(markers);
    UNUSED_VARIABLE(invGridSpacingSqr);
    
    vox::Size3 size = markers.size();
    
    bool hasX = false;
    double phiX = vox::kMaxD;
    
    if (i > 0) {
        if (vox::isInsideSdf(sign * output->tree().getValue(openvdb::Coord(i-1, j, k)) )) {
            hasX = true;
            phiX = std::min(phiX, sign * output->tree().getValue(openvdb::Coord(i-1, j, k)) );
        }
    }
    
    if (i + 1 < size.x) {
        if (vox::isInsideSdf(sign * output->tree().getValue(openvdb::Coord(i+1, j, k)) )) {
            hasX = true;
            phiX = std::min(phiX, sign * output->tree().getValue(openvdb::Coord(i+1, j, k)) );
        }
    }
    
    bool hasY = false;
    double phiY = vox::kMaxD;
    
    if (j > 0) {
        if (vox::isInsideSdf(sign * output->tree().getValue(openvdb::Coord(i, j-1, k)) )) {
            hasY = true;
            phiY = std::min(phiY, sign * output->tree().getValue(openvdb::Coord(i, j-1, k)) );
        }
    }
    
    if (j + 1 < size.y) {
        if (vox::isInsideSdf(sign * output->tree().getValue(openvdb::Coord(i, j+1, k)) )) {
            hasY = true;
            phiY = std::min(phiY, sign * output->tree().getValue(openvdb::Coord(i, j+1, k)) );
        }
    }
    
    bool hasZ = false;
    double phiZ = vox::kMaxD;
    
    if (k > 0) {
        if (vox::isInsideSdf(sign * output->tree().getValue(openvdb::Coord(i, j, k-1)) )) {
            hasZ = true;
            phiZ = std::min(phiZ, sign * output->tree().getValue(openvdb::Coord(i, j, k-1)) );
        }
    }
    
    if (k + 1 < size.z) {
        if (vox::isInsideSdf(sign * output->tree().getValue(openvdb::Coord(i, j, k+1)) )) {
            hasZ = true;
            phiZ = std::min(phiZ, sign * output->tree().getValue(openvdb::Coord(i, j, k+1)) );
        }
    }
    
    JET_ASSERT(hasX || hasY || hasZ);
    
    const double absCenter = std::fabs(output->tree().getValue(openvdb::Coord(i, j, k)) );
    
    double distToBndX =
    gridSpacing.x * absCenter / (absCenter + std::abs(phiX));
    
    double distToBndY =
    gridSpacing.y * absCenter / (absCenter + std::abs(phiY));
    
    double distToBndZ =
    gridSpacing.z * absCenter / (absCenter + std::abs(phiZ));
    
    double solution;
    double denomSqr = 0.0;
    
    if (hasX) {
        denomSqr += 1.0 / vox::square(distToBndX);
    }
    if (hasY) {
        denomSqr += 1.0 / vox::square(distToBndY);
    }
    if (hasZ) {
        denomSqr += 1.0 / vox::square(distToBndZ);
    }
    
    solution = 1.0 / std::sqrt(denomSqr);
    
    return sign * solution;
}

inline double solveQuad(const vox::Array3<char>& markers,
                        openvdb::DoubleGrid::Ptr output,
                        const vox::Vector3D& gridSpacing,
                        const vox::Vector3D& invGridSpacingSqr,
                        uint i, uint j, uint k) {
    vox::Size3 size = markers.size();
    
    bool hasX = false;
    double phiX = vox::kMaxD;
    
    if (i > 0) {
        if (markers(i - 1, j, k) == kKnown ) {
            hasX = true;
            phiX = std::min(phiX, output->tree().getValue(openvdb::Coord(i - 1, j, k)) );
        }
    }
    
    if (i + 1 < size.x) {
        if (markers(i + 1, j, k) == kKnown ) {
            hasX = true;
            phiX = std::min(phiX, output->tree().getValue(openvdb::Coord(i + 1, j, k)) );
        }
    }
    
    bool hasY = false;
    double phiY = vox::kMaxD;
    
    if (j > 0) {
        if (markers(i, j - 1, k) == kKnown ) {
            hasY = true;
            phiY = std::min(phiY, output->tree().getValue(openvdb::Coord(i, j - 1, k)) );
        }
    }
    
    if (j + 1 < size.y) {
        if (markers(i, j + 1, k) == kKnown ) {
            hasY = true;
            phiY = std::min(phiY, output->tree().getValue(openvdb::Coord(i, j + 1, k)) );
        }
    }
    
    bool hasZ = false;
    double phiZ = vox::kMaxD;
    
    if (k > 0) {
        if (markers(i, j, k - 1) == kKnown ) {
            hasZ = true;
            phiZ = std::min(phiZ, output->tree().getValue(openvdb::Coord(i, j, k - 1)) );
        }
    }
    
    if (k + 1 < size.z) {
        if (markers(i, j, k + 1) == kKnown ) {
            hasZ = true;
            phiZ = std::min(phiZ, output->tree().getValue(openvdb::Coord(i, j, k + 1)) );
        }
    }
    
    JET_ASSERT(hasX || hasY || hasZ);
    
    double solution = 0.0;
    
    // Initial guess
    if (hasX) {
        solution = std::max(solution, phiX + gridSpacing.x);
    }
    if (hasY) {
        solution = std::max(solution, phiY + gridSpacing.y);
    }
    if (hasZ) {
        solution = std::max(solution, phiZ + gridSpacing.z);
    }
    
    // Solve quad
    double a = 0.0;
    double b = 0.0;
    double c = -1.0;
    
    if (hasX) {
        a += invGridSpacingSqr.x;
        b -= phiX * invGridSpacingSqr.x;
        c += vox::square(phiX) * invGridSpacingSqr.x;
    }
    if (hasY) {
        a += invGridSpacingSqr.y;
        b -= phiY * invGridSpacingSqr.y;
        c += vox::square(phiY) * invGridSpacingSqr.y;
    }
    if (hasZ) {
        a += invGridSpacingSqr.z;
        b -= phiZ * invGridSpacingSqr.z;
        c += vox::square(phiZ) * invGridSpacingSqr.z;
    }
    
    double det = b * b - a * c;
    
    if (det > 0.0) {
        solution = (-b + std::sqrt(det)) / a;
    }
    
    return solution;
}

FmmLevelSetSolver3::FmmLevelSetSolver3() {}

void FmmLevelSetSolver3::reinitialize(const ScalarGrid3& inputSdf,
                                      double maxDistance,
                                      ScalarGrid3* outputSdf) {
    JET_THROW_INVALID_ARG_IF(!inputSdf.hasSameShape(*outputSdf));
    
    vox::Size3 size = inputSdf.dataSize();
    vox::Vector3D gridSpacing = inputSdf.gridSpacing();
    vox::Vector3D invGridSpacing = 1.0 / gridSpacing;
    vox::Vector3D invGridSpacingSqr = invGridSpacing * invGridSpacing;
    vox::Array3<char> markers(size);
    
    markers.forEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        outputSdf->getGrid()->tree().setValueOnly(coord, inputSdf.getGrid()->tree().getValue(coord));
    });
    
    auto output = outputSdf->getGrid()->treePtr();
    // Solve geometrically near the boundary
    markers.forEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(output->getValue(openvdb::Coord(i, j, k))) &&
            ((i > 0
              && !vox::isInsideSdf(output->getValue(openvdb::Coord(i - 1, j, k)))) ||
             (i + 1 < size.x
              && !vox::isInsideSdf(output->getValue(openvdb::Coord(i + 1, j, k)))) ||
             (j > 0
              && !vox::isInsideSdf(output->getValue(openvdb::Coord(i, j - 1, k)))) ||
             (j + 1 < size.y
              && !vox::isInsideSdf(output->getValue(openvdb::Coord(i, j + 1, k)))) ||
             (k > 0
              && !vox::isInsideSdf(output->getValue(openvdb::Coord(i, j, k - 1)))) ||
             (k + 1 < size.z
              && !vox::isInsideSdf(output->getValue(openvdb::Coord(i, j, k + 1)))) )) {
            output->setValueOnly(openvdb::Coord(i, j, k),
                                 solveQuadNearBoundary(
                                                       markers,
                                                       inputSdf.getGrid(),
                                                       gridSpacing,
                                                       invGridSpacingSqr,
                                                       -1.0, i, j, k));
        }
    });
    
    markers.forEachIndex([&](uint i, uint j, uint k) {
        if (!vox::isInsideSdf(output->getValue(openvdb::Coord(i, j, k))) &&
            ((i > 0
              && vox::isInsideSdf(output->getValue(openvdb::Coord(i - 1, j, k)))) ||
             (i + 1 < size.x
              && vox::isInsideSdf(output->getValue(openvdb::Coord(i + 1, j, k)))) ||
             (j > 0
              && vox::isInsideSdf(output->getValue(openvdb::Coord(i, j - 1, k)))) ||
             (j + 1 < size.y
              && vox::isInsideSdf(output->getValue(openvdb::Coord(i, j + 1, k)))) ||
             (k > 0
              && vox::isInsideSdf(output->getValue(openvdb::Coord(i, j, k - 1)))) ||
             (k + 1 < size.z
              && vox::isInsideSdf(output->getValue(openvdb::Coord(i, j, k + 1))))) ) {
            output->setValueOnly(openvdb::Coord(i, j, k),
                                 solveQuadNearBoundary(
                                                       markers,
                                                       inputSdf.getGrid(),
                                                       gridSpacing,
                                                       invGridSpacingSqr,
                                                       1.0, i, j, k));
        }
    });
    
    for (int sign = 0; sign < 2; ++sign) {
        // Build markers
        markers.parallelForEachIndex([&](uint i, uint j, uint k) {
            if (vox::isInsideSdf(output->getValue(openvdb::Coord(i, j, k)))) {
                markers(i, j, k) = kKnown;
            } else {
                markers(i, j, k) = kUnknown;
            }
        });
        
        auto compare = [&](const openvdb::Coord& a, const openvdb::Coord& b) {
            return output->getValue(a) > output->getValue(b);
        };
        
        // Enqueue initial candidates
        std::priority_queue<openvdb::Coord, std::vector<openvdb::Coord>, decltype(compare)>
        trial(compare);
        markers.forEachIndex([&](uint i, uint j, uint k) {
            if (markers(i, j, k) != kKnown &&
                ((i > 0 && markers(i - 1, j, k) == kKnown) ||
                 (i + 1 < size.x && markers(i + 1, j, k) == kKnown) ||
                 (j > 0 && markers(i, j - 1, k) == kKnown) ||
                 (j + 1 < size.y && markers(i, j + 1, k) == kKnown) ||
                 (k > 0 && markers(i, j, k - 1) == kKnown) ||
                 (k + 1 < size.z && markers(i, j, k + 1) == kKnown))) {
                trial.push(openvdb::Coord(i, j, k));
                markers(i, j, k) = kTrial;
            }
        });
        
        // Propagate
        while (!trial.empty()) {
            openvdb::Coord idx = trial.top();
            trial.pop();
            
            uint i = idx.x();
            uint j = idx.y();
            uint k = idx.z();
            
            markers(i, j, k) = kKnown;
            output->setValueOnly(idx, solveQuad(markers, outputSdf->getGrid(),
                                                gridSpacing,
                                                invGridSpacingSqr, i, j, k) );
            
            if (output->getValue(idx) > maxDistance) {
                break;
            }
            
            if (i > 0){
                if (markers(i - 1, j, k) == kUnknown) {
                    markers(i - 1, j, k) = kTrial;
                    output->setValueOnly(openvdb::Coord(i - 1, j, k),
                                         solveQuad(markers,
                                                   outputSdf->getGrid(),
                                                   gridSpacing,
                                                   invGridSpacingSqr,
                                                   i-1, j, k) );
                    trial.push(openvdb::Coord(i - 1, j, k));
                }
            }
            
            if (i + 1 < size.x) {
                if (markers(i + 1, j, k) == kUnknown) {
                    markers(i + 1, j, k) = kTrial;
                    output->setValueOnly(openvdb::Coord(i + 1, j, k),
                                         solveQuad(markers,
                                                   outputSdf->getGrid(),
                                                   gridSpacing,
                                                   invGridSpacingSqr,
                                                   i+1, j, k) );
                    trial.push(openvdb::Coord(i + 1, j, k));
                }
            }
            
            if (j > 0) {
                if (markers(i, j - 1, k) == kUnknown) {
                    markers(i, j - 1, k) = kTrial;
                    output->setValueOnly(openvdb::Coord(i, j-1, k),
                                         solveQuad(markers,
                                                   outputSdf->getGrid(),
                                                   gridSpacing,
                                                   invGridSpacingSqr,
                                                   i, j-1, k) );
                    trial.push(openvdb::Coord(i, j-1, k));
                }
            }
            
            if (j + 1 < size.y) {
                if (markers(i, j + 1, k) == kUnknown) {
                    markers(i, j + 1, k) = kTrial;
                    output->setValueOnly(openvdb::Coord(i, j+1, k),
                                         solveQuad(markers,
                                                   outputSdf->getGrid(),
                                                   gridSpacing,
                                                   invGridSpacingSqr,
                                                   i, j+1, k) );
                    trial.push(openvdb::Coord(i, j+1, k));
                }
            }
            
            if (k > 0) {
                if (markers(i, j, k - 1) == kUnknown) {
                    markers(i, j, k - 1) = kTrial;
                    output->setValueOnly(openvdb::Coord(i, j, k-1),
                                         solveQuad(markers,
                                                   outputSdf->getGrid(),
                                                   gridSpacing,
                                                   invGridSpacingSqr,
                                                   i, j, k-1) );
                    trial.push(openvdb::Coord(i, j, k-1));
                }
            }
            
            if (k + 1 < size.z) {
                if (markers(i, j, k + 1) == kUnknown) {
                    markers(i, j, k + 1) = kTrial;
                    output->setValueOnly(openvdb::Coord(i, j, k+1),
                                         solveQuad(markers,
                                                   outputSdf->getGrid(),
                                                   gridSpacing,
                                                   invGridSpacingSqr,
                                                   i, j, k+1) );
                    trial.push(openvdb::Coord(i, j, k+1));
                }
            }
        }
        
        // Flip the sign
        markers.parallelForEachIndex([&](uint i, uint j, uint k) {
            openvdb::Coord coord(i, j, k);
            output->setValueOnly(coord, -output->getValue(coord));
        });
    }
}

void FmmLevelSetSolver3::extrapolate(const ScalarGrid3& input,
                                     const vox::ScalarField3& sdf,
                                     double maxDistance, ScalarGrid3* output) {
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

void FmmLevelSetSolver3::extrapolate(const CollocatedVectorGrid3& input,
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

void FmmLevelSetSolver3::extrapolate(const FaceCenteredGrid3& input,
                                     const vox::ScalarField3& sdf,
                                     double maxDistance,
                                     FaceCenteredGrid3* output) {
    JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));
    
    const vox::Vector3D gridSpacing = input.gridSpacing();
    auto uPos = input.uPosition();
    
    vox::Array3<double> sdfAtU(input.uSize());
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

void FmmLevelSetSolver3::extrapolate(const openvdb::DoubleGrid::Ptr& input,
                                     const vox::ConstArrayAccessor3<double>& sdf,
                                     const vox::Vector3D& gridSpacing,
                                     double maxDistance,
                                     openvdb::DoubleGrid::Ptr output) {
    vox::Size3 size = sdf.size();
    vox::Vector3D invGridSpacing = 1.0 / gridSpacing;
    
    // Build markers
    vox::Array3<char> markers(size, kUnknown);
    markers.forEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(sdf(i, j, k))) {
            markers(i, j, k) = kKnown;
        }
        openvdb::Coord coord(i, j, k);
        output->tree().setValueOnly(coord, input->tree().getValue(coord) );
    });
    
    auto compare = [&](const openvdb::Coord& a, const openvdb::Coord& b) {
        return sdf(a.x(), a.y(), a.z()) > sdf(b.x(), b.y(), b.z());
    };
    
    // Enqueue initial candidates
    std::priority_queue<openvdb::Coord, std::vector<openvdb::Coord>, decltype(compare)>
    trial(compare);
    
    markers.forEachIndex([&](uint i, uint j, uint k) {
        if (markers(i, j, k) == kKnown) {
            return;
        }
        
        openvdb::Coord coord(i, j, k);
        if (i > 0 && markers(i - 1, j, k) == kKnown) {
            trial.push(coord);
            markers(i, j, k) = kTrial;
            return;
        }
        
        if (i + 1 < size.x && markers(i + 1, j, k) == kKnown) {
            trial.push(coord);
            markers(i, j, k) = kTrial;
            return;
        }
        
        if (j > 0 && markers(i, j - 1, k) == kKnown) {
            trial.push(coord);
            markers(i, j, k) = kTrial;
            return;
        }
        
        if (j + 1 < size.y && markers(i, j + 1, k) == kKnown) {
            trial.push(coord);
            markers(i, j, k) = kTrial;
            return;
        }
        
        if (k > 0 && markers(i, j, k - 1) == kKnown) {
            trial.push(coord);
            markers(i, j, k) = kTrial;
            return;
        }
        
        if (k + 1 < size.z && markers(i, j, k + 1) == kKnown) {
            trial.push(coord);
            markers(i, j, k) = kTrial;
            return;
        }
    });
    
    // Propagate
    while (!trial.empty()) {
        openvdb::Coord idx = trial.top();
        trial.pop();
        
        size_t i = idx.x();
        size_t j = idx.y();
        size_t k = idx.z();
        
        if (sdf(i, j, k) > maxDistance) {
            break;
        }
        
        vox::Vector3D grad = gradient3(sdf, gridSpacing,
                                       i, j, k).normalized();
        
        double sum = 0.0;
        double count = 0.0;
        
        openvdb::Coord neigh = idx + openvdb::Coord(-1, 0, 0);
        if (i > 0) {
            if (markers(i - 1, j, k) == kKnown) {
                double weight = std::max(grad.x, 0.0) * invGridSpacing.x;
                
                // If gradient is zero, then just assign 1 to weight
                if (weight < vox::kEpsilonD) {
                    weight = 1.0;
                }
                
                sum += weight * output->tree().getValue(neigh);
                count += weight;
            } else if (markers(i - 1, j, k) == kUnknown) {
                markers(i - 1, j, k) = kTrial;
                trial.push(neigh);
            }
        }
        
        neigh = idx + openvdb::Coord(1, 0, 0);
        if (i + 1 < size.x) {
            if (markers(i + 1, j, k) == kKnown) {
                double weight = std::max(grad.x, 0.0) * invGridSpacing.x;
                
                // If gradient is zero, then just assign 1 to weight
                if (weight < vox::kEpsilonD) {
                    weight = 1.0;
                }
                
                sum += weight * output->tree().getValue(neigh);
                count += weight;
            } else if (markers(i + 1, j, k) == kUnknown) {
                markers(i + 1, j, k) = kTrial;
                trial.push(neigh);
            }
        }
        
        neigh = idx + openvdb::Coord(0, -1, 0);
        if (j > 0) {
            if (markers(i, j - 1, k) == kKnown) {
                double weight = std::max(grad.y, 0.0) * invGridSpacing.y;
                
                // If gradient is zero, then just assign 1 to weight
                if (weight < vox::kEpsilonD) {
                    weight = 1.0;
                }
                
                sum += weight * output->tree().getValue(neigh);
                count += weight;
            } else if (markers(i, j - 1, k) == kUnknown) {
                markers(i, j - 1, k) = kTrial;
                trial.push(neigh);
            }
        }
        
        neigh = idx + openvdb::Coord(0, 1, 0);
        if (j + 1 < size.y) {
            if (markers(i, j + 1, k) == kKnown) {
                double weight = std::max(grad.y, 0.0) * invGridSpacing.y;
                
                // If gradient is zero, then just assign 1 to weight
                if (weight < vox::kEpsilonD) {
                    weight = 1.0;
                }
                
                sum += weight * output->tree().getValue(neigh);
                count += weight;
            } else if (markers(i, j + 1, k) == kUnknown) {
                markers(i, j + 1, k) = kTrial;
                trial.push(neigh);
            }
        }
        
        neigh = idx + openvdb::Coord(0, 0, -1);
        if (k > 0) {
            if (markers(i, j, k - 1) == kKnown) {
                double weight = std::max(grad.z, 0.0) * invGridSpacing.z;
                
                // If gradient is zero, then just assign 1 to weight
                if (weight < vox::kEpsilonD) {
                    weight = 1.0;
                }
                
                sum += weight * output->tree().getValue(neigh);
                count += weight;
            } else if (markers(i, j, k - 1) == kUnknown) {
                markers(i, j, k - 1) = kTrial;
                trial.push(neigh);
            }
        }
        
        neigh = idx + openvdb::Coord(0, 0, 1);
        if (k + 1 < size.z) {
            if (markers(i, j, k + 1) == kKnown) {
                double weight = std::max(grad.z, 0.0) * invGridSpacing.z;
                
                // If gradient is zero, then just assign 1 to weight
                if (weight < vox::kEpsilonD) {
                    weight = 1.0;
                }
                
                sum += weight * output->tree().getValue(neigh);
                count += weight;
            } else if (markers(i, j, k + 1) == kUnknown) {
                markers(i, j, k + 1) = kTrial;
                trial.push(neigh);
            }
        }
        
        JET_ASSERT(count > 0.0);
        
        output->tree().setValueOnly(idx, sum/count);
        markers(i, j, k) = kKnown;
    }
}
