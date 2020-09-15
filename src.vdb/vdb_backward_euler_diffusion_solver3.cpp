//
//  vdb_backward_euler_diffusion_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "vdb_backward_euler_diffusion_solver3.hpp"
#include "../src.common/constants.h"
#include "../src.common/fdm_iccg_solver3.h"
#include "../src.common/level_set_utils.h"
#include <openvdb/tools/PoissonSolver.h>

using namespace vdb;

const char kFluid = 0;
const char kAir = 1;
const char kBoundary = 2;

BackwardEulerDiffusionSolver3::BackwardEulerDiffusionSolver3(
                                                             BoundaryType boundaryType) :
_boundaryType(boundaryType) {
    _systemSolver = std::make_shared<vox::FdmIccgSolver3>(100, vox::kEpsilonD);
}

void BackwardEulerDiffusionSolver3::buildMarkers(
                                                 const vox::Size3& size,
                                                 const std::function<vox::Vector3D(const openvdb::Coord& coord)>& pos,
                                                 const vox::ScalarField3& boundarySdf,
                                                 const vox::ScalarField3& fluidSdf) {
    _markers.resize(size);
    
    _markers.parallelForEachIndex([&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(boundarySdf.sample(pos(openvdb::Coord(i, j, k))))) {
            _markers(i, j, k) = kBoundary;
        } else if (vox::isInsideSdf(fluidSdf.sample(pos(openvdb::Coord(i, j, k))))) {
            _markers(i, j, k) = kFluid;
        } else {
            _markers(i, j, k) = kAir;
        }
    });
}

void BackwardEulerDiffusionSolver3::buildMatrix(
                                                const vox::Size3& size,
                                                const vox::Vector3D& c) {
    _system.A.resize(size);
    
    bool isDirichlet = (_boundaryType == Dirichlet);
    
    // Build linear system
    _system.A.parallelForEachIndex(
                                   [&](size_t i, size_t j, size_t k) {
        auto& row = _system.A(i, j, k);
        
        // Initialize
        row.center = 1.0;
        row.right = row.up = row.front = 0.0;
        
        if (_markers(i, j, k) == kFluid) {
            if (i + 1 < size.x) {
                if ((isDirichlet && _markers(i + 1, j, k) != kAir)
                    || _markers(i + 1, j, k) == kFluid) {
                    row.center += c.x;
                }
                
                if (_markers(i + 1, j, k) == kFluid) {
                    row.right -=  c.x;
                }
            }
            
            if (i > 0
                && ((isDirichlet && _markers(i - 1, j, k) != kAir)
                    || _markers(i - 1, j, k) == kFluid)) {
                row.center += c.x;
            }
            
            if (j + 1 < size.y) {
                if ((isDirichlet && _markers(i, j + 1, k) != kAir)
                    || _markers(i, j + 1, k) == kFluid) {
                    row.center += c.y;
                }
                
                if (_markers(i, j + 1, k) == kFluid) {
                    row.up -=  c.y;
                }
            }
            
            if (j > 0
                && ((isDirichlet && _markers(i, j - 1, k) != kAir)
                    || _markers(i, j - 1, k) == kFluid)) {
                row.center += c.y;
            }
            
            if (k + 1 < size.z) {
                if ((isDirichlet && _markers(i, j, k + 1) != kAir)
                    || _markers(i, j, k + 1) == kFluid) {
                    row.center += c.z;
                }
                
                if (_markers(i, j, k + 1) == kFluid) {
                    row.front -=  c.z;
                }
            }
            
            if (k > 0
                && ((isDirichlet && _markers(i, j, k - 1) != kAir)
                    || _markers(i, j, k - 1) == kFluid)) {
                row.center += c.z;
            }
        }
    });
}

void BackwardEulerDiffusionSolver3::buildVectors(
                                                 const openvdb::DoubleGrid::Ptr f,
                                                 const vox::Size3 size,
                                                 const vox::Vector3D& c) {
    _system.x.resize(size, 0.0);
    _system.b.resize(size, 0.0);
    
    // Build linear system
    _system.x.parallelForEachIndex(
                                   [&](uint i, uint j, uint k) {
        _system.b(i, j, k) = _system.x(i, j, k) = f->tree().getValue(openvdb::Coord(i, j, k));
        
        if (_boundaryType == Dirichlet && _markers(i, j, k) == kFluid) {
            if (i + 1 < size.x && _markers(i + 1, j, k) == kBoundary) {
                _system.b(i, j, k) += c.x * f->tree().getValue(openvdb::Coord(i + 1, j, k));
            }
            
            if (i > 0 && _markers(i - 1, j, k) == kBoundary) {
                _system.b(i, j, k) += c.x * f->tree().getValue(openvdb::Coord(i - 1, j, k));
            }
            
            if (j + 1 < size.y && _markers(i, j + 1, k) == kBoundary) {
                _system.b(i, j, k) += c.y * f->tree().getValue(openvdb::Coord(i, j + 1, k));
            }
            
            if (j > 0 && _markers(i, j - 1, k) == kBoundary) {
                _system.b(i, j, k) += c.y * f->tree().getValue(openvdb::Coord(i, j - 1, k));
            }
            
            if (k + 1 < size.z && _markers(i, j, k + 1) == kBoundary) {
                _system.b(i, j, k) += c.z * f->tree().getValue(openvdb::Coord(i, j, k + 1));
            }
            
            if (k > 0 && _markers(i, j, k - 1) == kBoundary) {
                _system.b(i, j, k) += c.z * f->tree().getValue(openvdb::Coord(i, j, k - 1));
            }
        }
    });
}

void BackwardEulerDiffusionSolver3::buildVectors(
                                                 const openvdb::Vec3dGrid::Ptr f,
                                                 const vox::Size3 size,
                                                 const vox::Vector3D& c,
                                                 uint component) {
    _system.x.resize(size, 0.0);
    _system.b.resize(size, 0.0);
    
    // Build linear system
    _system.x.parallelForEachIndex(
                                   [&](uint i, uint j, uint k) {
        _system.b(i, j, k) = _system.x(i, j, k) = f->tree().getValue(openvdb::Coord(i, j, k))[component];
        
        if (_boundaryType == Dirichlet && _markers(i, j, k) == kFluid) {
            if (i + 1 < size.x && _markers(i + 1, j, k) == kBoundary) {
                _system.b(i, j, k) += c.x * f->tree().getValue(openvdb::Coord(i + 1, j, k))[component];
            }
            
            if (i > 0 && _markers(i - 1, j, k) == kBoundary) {
                _system.b(i, j, k) += c.x * f->tree().getValue(openvdb::Coord(i - 1, j, k))[component];
            }
            
            if (j + 1 < size.y && _markers(i, j + 1, k) == kBoundary) {
                _system.b(i, j, k) += c.y * f->tree().getValue(openvdb::Coord(i, j + 1, k))[component];
            }
            
            if (j > 0 && _markers(i, j - 1, k) == kBoundary) {
                _system.b(i, j, k) += c.y * f->tree().getValue(openvdb::Coord(i, j - 1, k))[component];
            }
            
            if (k + 1 < size.z && _markers(i, j, k + 1) == kBoundary) {
                _system.b(i, j, k) += c.z * f->tree().getValue(openvdb::Coord(i, j, k + 1))[component];
            }
            
            if (k > 0 && _markers(i, j, k - 1) == kBoundary) {
                _system.b(i, j, k) += c.z * f->tree().getValue(openvdb::Coord(i, j, k - 1))[component];
            }
        }
    });
}
//------------------------------------------------------------------------

void BackwardEulerDiffusionSolver3::solve(
                                          const ScalarGrid3& source,
                                          double diffusionCoefficient,
                                          double timeIntervalInSeconds,
                                          ScalarGrid3* dest,
                                          const vox::ScalarField3& boundarySdf,
                                          const vox::ScalarField3& fluidSdf) {
    auto pos = source.dataPosition();
    vox::Vector3D h = source.gridSpacing();
    vox::Vector3D c = timeIntervalInSeconds * diffusionCoefficient / (h * h);
    
    buildMarkers(source.dataSize(), pos, boundarySdf, fluidSdf);
    buildMatrix(source.dataSize(), c);
    buildVectors(source.getGrid(), source.dataSize(), c);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        _systemSolver->solve(&_system);
        
        // Assign the solution
        _markers.forEachIndex([&](uint i, uint j, uint k) {
            openvdb::Coord coord(i, j, k);
            (*dest).getGrid()->tree().setValueOnly(coord, _system.x(i, j, k));
        });
    }
}


void BackwardEulerDiffusionSolver3::solve(
                                          const CollocatedVectorGrid3& source,
                                          double diffusionCoefficient,
                                          double timeIntervalInSeconds,
                                          CollocatedVectorGrid3* dest,
                                          const vox::ScalarField3& boundarySdf,
                                          const vox::ScalarField3& fluidSdf) {
    auto pos = source.dataPosition();
    vox::Vector3D h = source.gridSpacing();
    vox::Vector3D c = timeIntervalInSeconds * diffusionCoefficient / (h * h);
    
    buildMarkers(source.dataSize(), pos, boundarySdf, fluidSdf);
    buildMatrix(source.dataSize(), c);
    
    vox::Array3<openvdb::Vec3d> result(source.dataSize());
    
    // u
    buildVectors(source.getGrid(), source.dataSize(), c, 0);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        _systemSolver->solve(&_system);
        
        // Assign the solution
        _markers.forEachIndex([&](uint i, uint j, uint k) {
            result(i, j, k).x() = _system.x(i, j, k);
        });
    }
    
    // v
    buildVectors(source.getGrid(), source.dataSize(), c, 1);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        _systemSolver->solve(&_system);
        
        // Assign the solution
        _markers.forEachIndex([&](uint i, uint j, uint k) {
            result(i, j, k).y() = _system.x(i, j, k);
        });
    }
    
    // w
    buildVectors(source.getGrid(), source.dataSize(), c, 2);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        _systemSolver->solve(&_system);
        
        // Assign the solution
        _markers.forEachIndex([&](uint i, uint j, uint k) {
            result(i, j, k).z() = _system.x(i, j, k);
        });
    }
    
    _markers.forEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        dest->getGrid()->tree().setValueOnly(coord, result(i, j, k));
    });
}

void BackwardEulerDiffusionSolver3::solve(
                                          const FaceCenteredGrid3& source,
                                          double diffusionCoefficient,
                                          double timeIntervalInSeconds,
                                          FaceCenteredGrid3* dest,
                                          const vox::ScalarField3& boundarySdf,
                                          const vox::ScalarField3& fluidSdf) {
    vox::Vector3D h = source.gridSpacing();
    vox::Vector3D c = timeIntervalInSeconds * diffusionCoefficient / (h * h);
    
    // u
    auto uPos = source.uPosition();
    buildMarkers(source.uSize(), uPos, boundarySdf, fluidSdf);
    buildMatrix(source.uSize(), c);
    buildVectors(source.getUGrid(), source.uSize(), c);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        _systemSolver->solve(&_system);
        
        // Assign the solution
        _markers.forEachIndex([&](uint i, uint j, uint k) {
            openvdb::Coord coord(i, j, k);
            dest->getUGrid()->tree().setValueOnly(coord, _system.x(i, j, k));
        });
    }
    
    // v
    auto vPos = source.vPosition();
    buildMarkers(source.vSize(), vPos, boundarySdf, fluidSdf);
    buildMatrix(source.vSize(), c);
    buildVectors(source.getVGrid(), source.vSize(), c);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        _systemSolver->solve(&_system);
        
        // Assign the solution
        _markers.forEachIndex([&](uint i, uint j, uint k) {
            openvdb::Coord coord(i, j, k);
            dest->getVGrid()->tree().setValueOnly(coord, _system.x(i, j, k));
        });
    }
    
    // w
    auto wPos = source.wPosition();
    buildMarkers(source.wSize(), wPos, boundarySdf, fluidSdf);
    buildMatrix(source.wSize(), c);
    buildVectors(source.getWGrid(), source.wSize(), c);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        _systemSolver->solve(&_system);
        
        // Assign the solution
        _markers.forEachIndex([&](uint i, uint j, uint k) {
            openvdb::Coord coord(i, j, k);
            dest->getWGrid()->tree().setValueOnly(coord, _system.x(i, j, k));
        });
    }
}
