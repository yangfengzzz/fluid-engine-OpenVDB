//
//  vdb_single_phase_pressure_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "../src.common/constants.h"
#include "../src.common/fdm_iccg_solver3.h"
#include "vdb_blocked_boundary_condition_solver3.hpp"
#include "vdb_single_phase_pressure_solver3.hpp"
#include "../src.common/level_set_utils.h"

using namespace vdb;

const char kFluid = 0;
const char kAir = 1;
const char kBoundary = 2;

const double kDefaultTolerance = 1e-6;

namespace {

void buildSingleSystem(vox::FdmMatrix3* A, vox::FdmVector3* b,
                       const vox::Array3<char>& markers,
                       const FaceCenteredGrid3& input) {
    vox::Size3 size = input.resolution();
    vox::Vector3D invH = 1.0 / input.gridSpacing();
    vox::Vector3D invHSqr = invH * invH;
    
    // Build linear system
    A->parallelForEachIndex([&](uint i, uint j, uint k) {
        auto& row = (*A)(i, j, k);
        
        // initialize
        row.center = row.right = row.up = row.front = 0.0;
        (*b)(i, j, k) = 0.0;
        
        if (markers(i, j, k) == kFluid) {
            (*b)(i, j, k) = input.divergenceAtCellCenter(i, j, k);
            
            if (i + 1 < size.x && markers(i + 1, j, k) != kBoundary) {
                row.center += invHSqr.x;
                if (markers(i + 1, j, k) == kFluid) {
                    row.right -= invHSqr.x;
                }
            }
            
            if (i > 0 && markers(i - 1, j, k) != kBoundary) {
                row.center += invHSqr.x;
            }
            
            if (j + 1 < size.y && markers(i, j + 1, k) != kBoundary) {
                row.center += invHSqr.y;
                if (markers(i, j + 1, k) == kFluid) {
                    row.up -= invHSqr.y;
                }
            }
            
            if (j > 0 && markers(i, j - 1, k) != kBoundary) {
                row.center += invHSqr.y;
            }
            
            if (k + 1 < size.z && markers(i, j, k + 1) != kBoundary) {
                row.center += invHSqr.z;
                if (markers(i, j, k + 1) == kFluid) {
                    row.front -= invHSqr.z;
                }
            }
            
            if (k > 0 && markers(i, j, k - 1) != kBoundary) {
                row.center += invHSqr.z;
            }
        } else {
            row.center = 1.0;
        }
    });
}

void buildSingleSystem(vox::MatrixCsrD* A, vox::VectorND* x, vox::VectorND* b,
                       const vox::Array3<char>& markers,
                       const FaceCenteredGrid3& input) {
    vox::Size3 size = input.resolution();
    vox::Vector3D invH = 1.0 / input.gridSpacing();
    vox::Vector3D invHSqr = invH * invH;
    
    const auto markerAcc = markers.constAccessor();
    
    A->clear();
    b->clear();
    
    size_t numRows = 0;
    vox::Array3<size_t> coordToIndex(size);
    markers.forEachIndex([&](size_t i, size_t j, size_t k) {
        const size_t cIdx = markerAcc.index(i, j, k);
        
        if (markerAcc[cIdx] == kFluid) {
            coordToIndex[cIdx] = numRows++;
        }
    });
    
    markers.forEachIndex([&](uint i, uint j, uint k) {
        const size_t cIdx = markerAcc.index(i, j, k);
        
        if (markerAcc[cIdx] == kFluid) {
            b->append(input.divergenceAtCellCenter(i, j, k));
            
            std::vector<double> row(1, 0.0);
            std::vector<size_t> colIdx(1, coordToIndex[cIdx]);
            
            if (i + 1 < size.x && markers(i + 1, j, k) != kBoundary) {
                row[0] += invHSqr.x;
                const size_t rIdx = markerAcc.index(i + 1, j, k);
                if (markers[rIdx] == kFluid) {
                    row.push_back(-invHSqr.x);
                    colIdx.push_back(coordToIndex[rIdx]);
                }
            }
            
            if (i > 0 && markers(i - 1, j, k) != kBoundary) {
                row[0] += invHSqr.x;
                const size_t lIdx = markerAcc.index(i - 1, j, k);
                if (markers[lIdx] == kFluid) {
                    row.push_back(-invHSqr.x);
                    colIdx.push_back(coordToIndex[lIdx]);
                }
            }
            
            if (j + 1 < size.y && markers(i, j + 1, k) != kBoundary) {
                row[0] += invHSqr.y;
                const size_t uIdx = markerAcc.index(i, j + 1, k);
                if (markers[uIdx] == kFluid) {
                    row.push_back(-invHSqr.y);
                    colIdx.push_back(coordToIndex[uIdx]);
                }
            }
            
            if (j > 0 && markers(i, j - 1, k) != kBoundary) {
                row[0] += invHSqr.y;
                const size_t dIdx = markerAcc.index(i, j - 1, k);
                if (markers[dIdx] == kFluid) {
                    row.push_back(-invHSqr.y);
                    colIdx.push_back(coordToIndex[dIdx]);
                }
            }
            
            if (k + 1 < size.z && markers(i, j, k + 1) != kBoundary) {
                row[0] += invHSqr.z;
                const size_t fIdx = markerAcc.index(i, j, k + 1);
                if (markers[fIdx] == kFluid) {
                    row.push_back(-invHSqr.z);
                    colIdx.push_back(coordToIndex[fIdx]);
                }
            }
            
            if (k > 0 && markers(i, j, k - 1) != kBoundary) {
                row[0] += invHSqr.z;
                const size_t bIdx = markerAcc.index(i, j, k - 1);
                if (markers[bIdx] == kFluid) {
                    row.push_back(-invHSqr.z);
                    colIdx.push_back(coordToIndex[bIdx]);
                }
            }
            
            A->addRow(row, colIdx);
        }
    });
    
    x->resize(b->size(), 0.0);
}

}

SinglePhasePressureSolver3::SinglePhasePressureSolver3() {
    _systemSolver = std::make_shared<vox::FdmIccgSolver3>(100, kDefaultTolerance);
}

SinglePhasePressureSolver3::~SinglePhasePressureSolver3() {}

BoundaryConditionSolver3Ptr
SinglePhasePressureSolver3::suggestedBoundaryConditionSolver() const {
    return std::make_shared<BlockedBoundaryConditionSolver3>();
}

const vox::FdmLinearSystemSolver3Ptr&
SinglePhasePressureSolver3::linearSystemSolver() const {
    return _systemSolver;
}

void SinglePhasePressureSolver3::setLinearSystemSolver(
                                                       const vox::FdmLinearSystemSolver3Ptr& solver) {
    _systemSolver = solver;
    _mgSystemSolver = std::dynamic_pointer_cast<vox::FdmMgSolver3>(_systemSolver);
    
    if (_mgSystemSolver == nullptr) {
        // In case of non-mg system, use flat structure.
        _mgSystem.clear();
    } else {
        // In case of mg system, use multi-level structure.
        _system.clear();
        _compSystem.clear();
    }
}

const vox::FdmVector3& SinglePhasePressureSolver3::pressure() const {
    if (_mgSystemSolver == nullptr) {
        return _system.x;
    } else {
        return _mgSystem.x.levels.front();
    }
}

void SinglePhasePressureSolver3::buildMarkers(const vox::Size3& size,
                                              const std::function<vox::Vector3D(const openvdb::Coord& coord)>& pos,
                                              const vox::ScalarField3& boundarySdf,
                                              const vox::ScalarField3& fluidSdf){
    // Build levels
    size_t maxLevels = 1;
    if (_mgSystemSolver != nullptr) {
        maxLevels = _mgSystemSolver->params().maxNumberOfLevels;
    }
    vox::FdmMgUtils3::resizeArrayWithFinest(size, maxLevels, &_markers);
    
    // Build top-level markers
    _markers[0].parallelForEachIndex([&](uint i, uint j, uint k) {
        vox::Vector3D pt = pos(openvdb::Coord(i, j, k));
        if (vox::isInsideSdf(boundarySdf.sample(pt))) {
            _markers[0](i, j, k) = kBoundary;
        } else if (vox::isInsideSdf(fluidSdf.sample(pt))) {
            _markers[0](i, j, k) = kFluid;
        } else {
            _markers[0](i, j, k) = kAir;
        }
    });
    
    // Build sub-level markers
    for (size_t l = 1; l < _markers.size(); ++l) {
        const auto& finer = _markers[l - 1];
        auto& coarser = _markers[l];
        const vox::Size3 n = coarser.size();
        
        vox::parallelRangeFor(
                              vox::kZeroSize, n.x,
                              vox::kZeroSize, n.y,
                              vox::kZeroSize, n.z,
                              [&](size_t iBegin, size_t iEnd, size_t jBegin, size_t jEnd,
                                  size_t kBegin, size_t kEnd) {
            std::array<size_t, 4> kIndices;
            
            for (size_t k = kBegin; k < kEnd; ++k) {
                kIndices[0] = (k > 0) ? 2 * k - 1 : 2 * k;
                kIndices[1] = 2 * k;
                kIndices[2] = 2 * k + 1;
                kIndices[3] = (k + 1 < n.z) ? 2 * k + 2 : 2 * k + 1;
                
                std::array<size_t, 4> jIndices;
                
                for (size_t j = jBegin; j < jEnd; ++j) {
                    jIndices[0] = (j > 0) ? 2 * j - 1 : 2 * j;
                    jIndices[1] = 2 * j;
                    jIndices[2] = 2 * j + 1;
                    jIndices[3] = (j + 1 < n.y) ? 2 * j + 2 : 2 * j + 1;
                    
                    std::array<size_t, 4> iIndices;
                    for (size_t i = iBegin; i < iEnd; ++i) {
                        iIndices[0] = (i > 0) ? 2 * i - 1 : 2 * i;
                        iIndices[1] = 2 * i;
                        iIndices[2] = 2 * i + 1;
                        iIndices[3] = (i + 1 < n.x) ? 2 * i + 2 : 2 * i + 1;
                        
                        int cnt[3] = {0, 0, 0};
                        for (size_t z = 0; z < 4; ++z) {
                            for (size_t y = 0; y < 4; ++y) {
                                for (size_t x = 0; x < 4; ++x) {
                                    char f = finer(iIndices[x], jIndices[y],
                                                   kIndices[z]);
                                    if (f == kBoundary) {
                                        ++cnt[(int)kBoundary];
                                    } else if (f == kFluid) {
                                        ++cnt[(int)kFluid];
                                    } else {
                                        ++cnt[(int)kAir];
                                    }
                                }
                            }
                        }
                        
                        coarser(i, j, k) = static_cast<char>(
                                                             vox::argmax3(cnt[0], cnt[1], cnt[2]));
                    }
                }
            }
        });
    }
}

void SinglePhasePressureSolver3::decompressSolution() {
    const auto acc = _markers[0].constAccessor();
    _system.x.resize(acc.size());
    
    size_t row = 0;
    _markers[0].forEachIndex([&](size_t i, size_t j, size_t k) {
        if (acc(i, j, k) == kFluid) {
            _system.x(i, j, k) = _compSystem.x[row];
            ++row;
        }
    });
}

void SinglePhasePressureSolver3::buildSystem(const FaceCenteredGrid3& input,
                                             bool useCompressed) {
    vox::Size3 size = input.resolution();
    size_t numLevels = 1;
    
    if (_mgSystemSolver == nullptr) {
        if (!useCompressed) {
            _system.resize(size);
        }
    } else {
        // Build levels
        size_t maxLevels = _mgSystemSolver->params().maxNumberOfLevels;
        vox::FdmMgUtils3::resizeArrayWithFinest(size, maxLevels,
                                                &_mgSystem.A.levels);
        vox::FdmMgUtils3::resizeArrayWithFinest(size, maxLevels,
                                                &_mgSystem.x.levels);
        vox::FdmMgUtils3::resizeArrayWithFinest(size, maxLevels,
                                                &_mgSystem.b.levels);
        
        numLevels = _mgSystem.A.levels.size();
    }
    
    // Build top level
    const FaceCenteredGrid3* finer = &input;
    if (_mgSystemSolver == nullptr) {
        if (useCompressed) {
            buildSingleSystem(&_compSystem.A, &_compSystem.x, &_compSystem.b,
                              _markers[0], *finer);
        } else {
            buildSingleSystem(&_system.A, &_system.b, _markers[0], *finer);
        }
    } else {
        buildSingleSystem(&_mgSystem.A.levels.front(),
                          &_mgSystem.b.levels.front(), _markers[0], *finer);
    }
    
    // Build sub-levels
    FaceCenteredGrid3 coarser;
    for (size_t l = 1; l < numLevels; ++l) {
        auto res = finer->resolution();
        auto h = finer->gridSpacing();
        auto o = finer->origin();
        res.x = res.x >> 1;
        res.y = res.y >> 1;
        res.z = res.z >> 1;
        h *= 2.0;
        
        // Down sample
        coarser.resize(res, h, o);
        coarser.fill(finer->sampler());
        
        buildSingleSystem(&_mgSystem.A.levels[l], &_mgSystem.b.levels[l],
                          _markers[l], coarser);
        
        finer = &coarser;
    }
}

void SinglePhasePressureSolver3::applyPressureGradient(
                                                       const FaceCenteredGrid3& input,
                                                       FaceCenteredGrid3* output) {
    vox::Size3 size = input.resolution();
    auto uPos = input.uPosition();
    auto vPos = input.vPosition();
    auto wPos = input.wPosition();
    auto u = input.getUGrid()->tree();
    auto v = input.getVGrid()->tree();
    auto w = input.getWGrid()->tree();
    
    const auto& x = pressure();
    
    vox::Vector3D invH = 1.0 / input.gridSpacing();
    
    x.forEachIndex([&](uint i, uint j, uint k) {
        if (_markers[0](i, j, k) == kFluid) {
            if (i + 1 < size.x && _markers[0](i + 1, j, k) != kBoundary) {
                output->getUGrid()->tree().setValue(openvdb::Coord(i + 1, j, k),
                                                    input.u(openvdb::Coord(i + 1, j, k))
                                                    + invH.x * (x(i + 1, j, k) - x(i, j, k)));
            }
            if (j + 1 < size.y && _markers[0](i, j + 1, k) != kBoundary) {
                output->getVGrid()->tree().setValue(openvdb::Coord(i, j + 1, k),
                                                    input.v(openvdb::Coord(i, j + 1, k))
                                                    + invH.y * (x(i, j + 1, k) - x(i, j, k)));
            }
            if (k + 1 < size.z && _markers[0](i, j, k + 1) != kBoundary) {
                output->getWGrid()->tree().setValue(openvdb::Coord(i, j, k + 1),
                                                    input.w(openvdb::Coord(i, j, k + 1))
                                                    + invH.z * (x(i, j, k + 1) - x(i, j, k)));
            }
        }
    });
}

void SinglePhasePressureSolver3::solve(const FaceCenteredGrid3& input,
                                       double timeIntervalInSeconds,
                                       FaceCenteredGrid3* output,
                                       const vox::ScalarField3& boundarySdf,
                                       const vox::VectorField3& boundaryVelocity,
                                       const vox::ScalarField3& fluidSdf,
                                       bool useCompressed) {
    UNUSED_VARIABLE(timeIntervalInSeconds);
    UNUSED_VARIABLE(boundaryVelocity);
    
    auto pos = input.cellCenterPosition();
    buildMarkers(input.resolution(), pos, boundarySdf, fluidSdf);
    buildSystem(input, useCompressed);
    
    if (_systemSolver != nullptr) {
        // Solve the system
        if (_mgSystemSolver == nullptr) {
            if (useCompressed) {
                _system.clear();
                _systemSolver->solveCompressed(&_compSystem);
                decompressSolution();
            } else {
                _compSystem.clear();
                _systemSolver->solve(&_system);
            }
        } else {
            _mgSystemSolver->solve(&_mgSystem);
        }
        
        // Apply pressure gradient
        applyPressureGradient(input, output);
    }
}
