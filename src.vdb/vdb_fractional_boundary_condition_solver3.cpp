//
//  vdb_fractional_boundary_condition_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_helper.h"
#include "../src.common/physics_helpers.h"
#include "../src.common/array_utils.h"
#include "vdb_fractional_boundary_condition_solver3.hpp"
#include "../src.common/level_set_utils.h"
#include "../src.common/surface_to_implicit3.h"
#include <openvdb/math/Stencils.h>
#include <algorithm>

using namespace vdb;

FractionalBoundaryConditionSolver3
::FractionalBoundaryConditionSolver3() {
}

FractionalBoundaryConditionSolver3::
~FractionalBoundaryConditionSolver3() {
}

vox::ScalarField3Ptr FractionalBoundaryConditionSolver3::colliderSdf() const {
    return _colliderSdf;
}

vox::VectorField3Ptr
FractionalBoundaryConditionSolver3::colliderVelocityField() const {
    return _colliderVel;
}

void FractionalBoundaryConditionSolver3::onColliderUpdated(const vox::Size3& gridSize,
                                                           const vox::Vector3D& gridSpacing,
                                                           const vox::Vector3D& gridOrigin) {
    if (_colliderSdf == nullptr) {
        _colliderSdf = std::make_shared<CellCenteredScalarGrid3>();
    }
    
    _colliderSdf->resize(gridSize, gridSpacing, gridOrigin);
    
    if (collider() != nullptr) {
        vox::Surface3Ptr surface = collider()->surface();
        vox::ImplicitSurface3Ptr implicitSurface
        = std::dynamic_pointer_cast<vox::ImplicitSurface3>(surface);
        if (implicitSurface == nullptr) {
            implicitSurface = std::make_shared<vox::SurfaceToImplicit3>(surface);
        }
        
        _colliderSdf->fill([&](const vox::Vector3D& pt) {
            return implicitSurface->signedDistance(pt);
        }, vox::ExecutionPolicy::kSerial);
        
        _colliderVel = vox::CustomVectorField3::builder()
        .withFunction([&] (const vox::Vector3D& x) {
            return collider()->velocityAt(x);
        })
        .withDerivativeResolution(gridSpacing.x)
        .makeShared();
    } else {
        _colliderSdf->fill(vox::kMaxD, vox::ExecutionPolicy::kSerial);
        
        _colliderVel = vox::CustomVectorField3::builder()
        .withFunction([] (const vox::Vector3D&) {
            return vox::Vector3D();
        })
        .withDerivativeResolution(gridSpacing.x)
        .makeShared();
    }
}

void FractionalBoundaryConditionSolver3::constrainVelocity(
                                                           FaceCenteredGrid3* velocity,
                                                           unsigned int extrapolationDepth) {
    vox::Size3 size = velocity->resolution();
    if (_colliderSdf == nullptr || _colliderSdf->resolution() != size) {
        updateCollider(
                       collider(),
                       size,
                       velocity->gridSpacing(),
                       velocity->origin());
    }
    
    auto uPos = velocity->uPosition();
    auto vPos = velocity->vPosition();
    auto wPos = velocity->wPosition();
    
    vox::Array3<char> uMarker(velocity->uSize(), 1);
    vox::Array3<char> vMarker(velocity->vSize(), 1);
    vox::Array3<char> wMarker(velocity->wSize(), 1);
    
    vox::Vector3D h = velocity->gridSpacing();
    
    // Assign collider's velocity first and initialize markers
    velocity->parallelForEachUIndex([&](const openvdb::Coord& coord) {
        vox::Vector3D pt = uPos(coord);
        double phi0 = _colliderSdf->sample(pt - vox::Vector3D(0.5 * h.x, 0.0, 0.0));
        double phi1 = _colliderSdf->sample(pt + vox::Vector3D(0.5 * h.x, 0.0, 0.0));
        double frac = vox::fractionInsideSdf(phi0, phi1);
        frac = 1.0 - vox::clamp(frac, 0.0, 1.0);
        
        if (frac > 0.0) {
            uMarker(coord.x(), coord.y(), coord.z()) = 1;
        } else {
            vox::Vector3D colliderVel = collider()->velocityAt(pt);
            velocity->getUGrid()->tree().setValue(coord, colliderVel.x);
            uMarker(coord.x(), coord.y(), coord.z()) = 0;
        }
    });
    
    velocity->parallelForEachVIndex([&](const openvdb::Coord& coord) {
        vox::Vector3D pt = vPos(coord);
        double phi0 = _colliderSdf->sample(pt - vox::Vector3D(0.0, 0.5 * h.y, 0.0));
        double phi1 = _colliderSdf->sample(pt + vox::Vector3D(0.0, 0.5 * h.y, 0.0));
        double frac = vox::fractionInsideSdf(phi0, phi1);
        frac = 1.0 - vox::clamp(frac, 0.0, 1.0);
        
        if (frac > 0.0) {
            vMarker(coord.x(), coord.y(), coord.z()) = 1;
        } else {
            vox::Vector3D colliderVel = collider()->velocityAt(pt);
            velocity->getVGrid()->tree().setValue(coord, colliderVel.y);
            vMarker(coord.x(), coord.y(), coord.z()) = 0;
        }
    });
    
    velocity->parallelForEachWIndex([&](const openvdb::Coord& coord) {
        vox::Vector3D pt = wPos(coord);
        double phi0 = _colliderSdf->sample(pt - vox::Vector3D(0.0, 0.0, 0.5 * h.z ));
        double phi1 = _colliderSdf->sample(pt + vox::Vector3D(0.0, 0.0, 0.5 * h.z ));
        double frac = vox::fractionInsideSdf(phi0, phi1);
        frac = 1.0 - vox::clamp(frac, 0.0, 1.0);
        
        if (frac > 0.0) {
            wMarker(coord.x(), coord.y(), coord.z()) = 1;
        } else {
            vox::Vector3D colliderVel = collider()->velocityAt(pt);
            velocity->getWGrid()->tree().setValue(coord, colliderVel.z);
            wMarker(coord.x(), coord.y(), coord.z()) = 0;
        }
    });
    
    // Free-slip: Extrapolate fluid velocity into the collider
    extrapolateToRegion<openvdb::DoubleGrid>(velocity->getUGrid(),
                                             uMarker,
                                             extrapolationDepth,
                                             velocity->getUGrid());
    extrapolateToRegion<openvdb::DoubleGrid>(velocity->getVGrid(),
                                             vMarker,
                                             extrapolationDepth,
                                             velocity->getVGrid());
    extrapolateToRegion<openvdb::DoubleGrid>(velocity->getWGrid(),
                                             wMarker,
                                             extrapolationDepth,
                                             velocity->getWGrid());
    
    // No-flux: project the extrapolated velocity to the collider's surface
    // normal
    openvdb::DoubleGrid::Ptr uTemp = openvdb::DoubleGrid::create(0.0);
    uTemp->topologyUnion(*velocity->getUGrid());
    velocity->parallelForEachUIndex([&](const openvdb::Coord& coord) {
        vox::Vector3D pt = uPos(coord);
        if (vox::isInsideSdf(_colliderSdf->sample(pt))) {
            vox::Vector3D colliderVel = collider()->velocityAt(pt);
            vox::Vector3D vel = velocity->sample(pt);
            vox::Vector3D g = _colliderSdf->gradient(pt);
            
            if (g.lengthSquared() > 0.0) {
                vox::Vector3D n = g.normalized();
                vox::Vector3D velr = vel - colliderVel;
                vox::Vector3D velt =
                vox::projectAndApplyFriction(
                                             velr, n,
                                             collider()->frictionCoefficient());
                
                vox::Vector3D velp = velt + colliderVel;
                uTemp->tree().setValueOn(coord, velp.x );
            } else {
                uTemp->tree().setValueOn(coord, colliderVel.x );
            }
        } else {
            uTemp->tree().setValueOn(coord,
                                     velocity->getUGrid()->tree().getValue(coord));
        }
    });
    
    openvdb::DoubleGrid::Ptr vTemp = openvdb::DoubleGrid::create(0.0);
    vTemp->topologyUnion(*velocity->getVGrid());
    velocity->parallelForEachVIndex([&](const openvdb::Coord& coord) {
        vox::Vector3D pt = vPos(coord);
        if (vox::isInsideSdf(_colliderSdf->sample(pt))) {
            vox::Vector3D colliderVel = collider()->velocityAt(pt);
            vox::Vector3D vel = velocity->sample(pt);
            vox::Vector3D g = _colliderSdf->gradient(pt);
            
            if (g.lengthSquared() > 0.0) {
                vox::Vector3D n = g.normalized();
                vox::Vector3D velr = vel - colliderVel;
                vox::Vector3D velt =
                vox::projectAndApplyFriction(
                                             velr, n,
                                             collider()->frictionCoefficient());
                
                vox::Vector3D velp = velt + colliderVel;
                vTemp->tree().setValueOn(coord, velp.y );
            } else {
                vTemp->tree().setValueOn(coord, colliderVel.y );
            }
        } else {
            vTemp->tree().setValueOn(coord,
                                     velocity->getVGrid()-> tree().getValue(coord));
        }
    });
    
    openvdb::DoubleGrid::Ptr wTemp = openvdb::DoubleGrid::create(0.0);
    wTemp->topologyUnion(*velocity->getWGrid());
    velocity->parallelForEachWIndex([&](const openvdb::Coord& coord) {
        vox::Vector3D pt = wPos(coord);
        if (vox::isInsideSdf(_colliderSdf->sample(pt))) {
            vox::Vector3D colliderVel = collider()->velocityAt(pt);
            vox::Vector3D vel = velocity->sample(pt);
            vox::Vector3D g = _colliderSdf->gradient(pt);
            
            if (g.lengthSquared() > 0.0) {
                vox::Vector3D n = g.normalized();
                vox::Vector3D velr = vel - colliderVel;
                vox::Vector3D velt =
                vox::projectAndApplyFriction(
                                             velr, n,
                                             collider()->frictionCoefficient());
                
                vox::Vector3D velp = velt + colliderVel;
                wTemp->tree().setValueOn(coord, velp.z );
            } else {
                wTemp->tree().setValueOn(coord, colliderVel.z );
            }
        } else {
            wTemp->tree().setValueOn(coord,
                                     velocity->getWGrid()-> tree().getValue(coord));
        }
    });
    
    // Transfer results
    velocity->parallelForEachUIndex([&](const openvdb::Coord& coord) {
        velocity->getUGrid()->tree().setValueOn(coord, uTemp->tree().getValue(coord));
    });
    velocity->parallelForEachVIndex([&](const openvdb::Coord& coord) {
        velocity->getVGrid()->tree().setValueOn(coord, vTemp->tree().getValue(coord));
    });
    velocity->parallelForEachWIndex([&](const openvdb::Coord& coord) {
        velocity->getWGrid()->tree().setValueOn(coord, wTemp->tree().getValue(coord));
    });
    
    if (closedDomainBoundaryFlag() & vox::kDirectionLeft) {
        for (uint k = 0; k < velocity->uSize().z; ++k) {
            for (uint j = 0; j < velocity->uSize().y; ++j) {
                velocity->getUGrid()->tree().setValueOnly(openvdb::Coord(0, j, k), 0);
            }
        }
    }
    if (closedDomainBoundaryFlag() & vox::kDirectionRight) {
        for (uint k = 0; k < velocity->uSize().z; ++k) {
            for (uint j = 0; j < velocity->uSize().y; ++j) {
                velocity->getUGrid()->tree().setValueOnly(openvdb::Coord((int)velocity->uSize().x - 1,
                                                                         j, k),
                                                          0);
            }
        }
    }
    if (closedDomainBoundaryFlag() & vox::kDirectionDown) {
        for (uint k = 0; k < velocity->vSize().z; ++k) {
            for (uint i = 0; i < velocity->vSize().x; ++i) {
                velocity->getVGrid()->tree().setValueOnly(openvdb::Coord(i,
                                                                         0, k),
                                                          0);
            }
        }
    }
    if (closedDomainBoundaryFlag() & vox::kDirectionUp) {
        for (uint k = 0; k < velocity->vSize().z; ++k) {
            for (uint i = 0; i < velocity->vSize().x; ++i) {
                velocity->getVGrid()->tree().setValueOnly(openvdb::Coord(i,
                                                                         (int)velocity->vSize().y - 1,
                                                                         k),
                                                          0);
            }
        }
    }
    if (closedDomainBoundaryFlag() & vox::kDirectionBack) {
        for (uint j = 0; j < velocity->wSize().y; ++j) {
            for (uint i = 0; i < velocity->wSize().x; ++i) {
                velocity->getWGrid()->tree().setValueOnly(openvdb::Coord(i,
                                                                         j, 0),
                                                          0);
            }
        }
    }
    if (closedDomainBoundaryFlag() & vox::kDirectionFront) {
        for (uint j = 0; j < velocity->wSize().y; ++j) {
            for (uint i = 0; i < velocity->wSize().x; ++i) {
                velocity->getWGrid()->tree().setValueOnly(openvdb::Coord(i, j,
                                                                         (int)velocity->wSize().z - 1),
                                                          0);
            }
        }
    }
}
