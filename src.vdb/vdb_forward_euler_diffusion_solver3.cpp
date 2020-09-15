//
//  vdb_forward_euler_diffusion_solver3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/7.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "../src.common/fdm_utils.h"
#include "vdb_fdm_utils.h"
#include "vdb_helper.h"
#include "vdb_forward_euler_diffusion_solver3.hpp"
#include "../src.common/level_set_utils.h"

using namespace vdb;

static const char kFluid = 0;
static const char kAir = 1;
static const char kBoundary = 2;

template <typename GridType>
typename GridType::ValueType laplacian(
                                       const typename GridType::Ptr& data,
                                       const vox::Array3<char>& marker,
                                       const vox::Vector3D& gridSpacing,
                                       openvdb::Coord ijk) {
    const typename GridType::ValueType center = data->tree().getValue(ijk);
    const vox::Size3 ds = marker.size();
    JET_ASSERT(ijk.x() < ds.x && ijk.y() < ds.y && ijk.z() < ds.z);
    
    typename GridType::ValueType dleft = openvdb::zeroVal<typename GridType::ValueType>();
    typename GridType::ValueType dright = openvdb::zeroVal<typename GridType::ValueType>();
    typename GridType::ValueType ddown = openvdb::zeroVal<typename GridType::ValueType>();
    typename GridType::ValueType dup = openvdb::zeroVal<typename GridType::ValueType>();
    typename GridType::ValueType dback = openvdb::zeroVal<typename GridType::ValueType>();
    typename GridType::ValueType dfront = openvdb::zeroVal<typename GridType::ValueType>();
    
    openvdb::Coord neigh = ijk + openvdb::Coord(-1, 0, 0);
    if (ijk.x() > 0
        && marker(neigh.x(), neigh.y(), neigh.z()) == kFluid) {
        dleft = center - data->tree().getValue(neigh);
    }
    
    neigh = ijk + openvdb::Coord(1, 0, 0);
    if (ijk.x() + 1 < ds.x
        && marker(neigh.x(), neigh.y(), neigh.z()) == kFluid) {
        dright = data->tree().getValue(neigh) - center;
    }
    
    neigh = ijk + openvdb::Coord(0, -1, 0);
    if (ijk.y() > 0
        && marker(neigh.x(), neigh.y(), neigh.z()) == kFluid) {
        ddown = center - data->tree().getValue(neigh);
    }
    
    neigh = ijk + openvdb::Coord(0, 1, 0);
    if (ijk.y() + 1 < ds.y
        && marker(neigh.x(), neigh.y(), neigh.z()) == kFluid) {
        dup = data->tree().getValue(neigh) - center;
    }
    
    neigh = ijk + openvdb::Coord(0, 0, -1);
    if (ijk.z() > 0
        && marker(neigh.x(), neigh.y(), neigh.z()) == kFluid) {
        dback = center - data->tree().getValue(neigh);
    }
    
    neigh = ijk + openvdb::Coord(0, 0, 1);
    if (ijk.z() + 1 < ds.z
        && marker(neigh.x(), neigh.y(), neigh.z()) == kFluid) {
        dfront = data->tree().getValue(neigh) - center;
    }
    
    return (dright - dleft) / vox::square(gridSpacing.x)
    + (dup - ddown) / vox::square(gridSpacing.y)
    + (dfront - dback) / vox::square(gridSpacing.z);
}

ForwardEulerDiffusionSolver3::ForwardEulerDiffusionSolver3() {
}

void ForwardEulerDiffusionSolver3::solve(
                                         const ScalarGrid3& source,
                                         double diffusionCoefficient,
                                         double timeIntervalInSeconds,
                                         ScalarGrid3* dest,
                                         const vox::ScalarField3& boundarySdf,
                                         const vox::ScalarField3& fluidSdf) {
    vox::Vector3D h = source.gridSpacing();
    auto pos = source.dataPosition();
    
    buildMarkers(source.resolution(), pos, boundarySdf, fluidSdf);
    
    _markers.forEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        if (_markers(i, j, k) == kFluid) {
            dest->getGrid()->tree().setValueOnly(coord,
                                                 source.getGrid()->tree().getValue(coord)
                                                 + diffusionCoefficient
                                                 * timeIntervalInSeconds
                                                 * laplacian<openvdb::DoubleGrid>(source.getGrid(),
                                                                                  _markers, h, coord) );
        } else {
            dest->getGrid()->tree().setValueOnly(coord,
                                                 source.getGrid()->tree().getValue(coord) );
        }
    });
}


void ForwardEulerDiffusionSolver3::solve(
                                         const CollocatedVectorGrid3& source,
                                         double diffusionCoefficient,
                                         double timeIntervalInSeconds,
                                         CollocatedVectorGrid3* dest,
                                         const vox::ScalarField3& boundarySdf,
                                         const vox::ScalarField3& fluidSdf) {
    vox::Vector3D h = source.gridSpacing();
    auto pos = source.dataPosition();
    
    buildMarkers(source.resolution(), pos, boundarySdf, fluidSdf);
    
    _markers.forEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        if (_markers(i, j, k) == kFluid) {
            dest->getGrid()->tree().setValueOnly(coord, source.getGrid()->tree().getValue(coord)
                                                 + diffusionCoefficient
                                                 * timeIntervalInSeconds
                                                 * laplacian<openvdb::Vec3DGrid>(source.getGrid(),
                                                                                 _markers, h, coord) );
        } else {
            dest->getGrid()->tree().setValueOnly(coord, source.getGrid()->tree().getValue(coord) );
        }
    });
}

void ForwardEulerDiffusionSolver3::solve(
                                         const FaceCenteredGrid3& source,
                                         double diffusionCoefficient,
                                         double timeIntervalInSeconds,
                                         FaceCenteredGrid3* dest,
                                         const vox::ScalarField3& boundarySdf,
                                         const vox::ScalarField3& fluidSdf) {
    auto uPos = source.uPosition();
    auto vPos = source.vPosition();
    auto wPos = source.wPosition();
    vox::Vector3D h = source.gridSpacing();
    
    buildMarkers(source.uSize(), uPos, boundarySdf, fluidSdf);
    
    _markers.forEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        if (!vox::isInsideSdf(boundarySdf.sample(source.uPosition()(coord) ))) {
            dest->getUGrid()->tree().setValueOnly(coord,
                                                  source.getUGrid()->tree().getValue(coord)
                                                  + diffusionCoefficient
                                                  * timeIntervalInSeconds
                                                  * laplacian3(source.getUGrid(),
                                                               source.uSize(), h,
                                                               coord.x(), coord.y(), coord.z()) );
        }
    });
    
    buildMarkers(source.vSize(), vPos, boundarySdf, fluidSdf);
    
    _markers.forEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        if (!vox::isInsideSdf(boundarySdf.sample(source.vPosition()(coord) ))) {
            dest->getVGrid()->tree().setValueOnly(coord,
                                                  source.getVGrid()->tree().getValue(coord)
                                                  + diffusionCoefficient
                                                  * timeIntervalInSeconds
                                                  * laplacian3(source.getVGrid(),
                                                               source.vSize(), h,
                                                               coord.x(), coord.y(), coord.z()) );
        }
    });
    
    buildMarkers(source.wSize(), wPos, boundarySdf, fluidSdf);
    
    _markers.forEachIndex([&](uint i, uint j, uint k) {
        openvdb::Coord coord(i, j, k);
        if (!vox::isInsideSdf(boundarySdf.sample(source.wPosition()(coord) ))) {
            dest->getWGrid()->tree().setValueOnly(coord,
                                                  source.getWGrid()->tree().getValue(coord)
                                                  + diffusionCoefficient
                                                  * timeIntervalInSeconds
                                                  * laplacian3(source.getWGrid(),
                                                               source.wSize(), h,
                                                               coord.x(), coord.y(), coord.z()) );
        }
    });
}

void ForwardEulerDiffusionSolver3::buildMarkers(
                                                const vox::Size3& size,
                                                const std::function<vox::Vector3D(const openvdb::Coord& coord)>& pos,
                                                const vox::ScalarField3& boundarySdf,
                                                const vox::ScalarField3& fluidSdf){
    _markers.resize(size);
    
    _markers.forEachIndex(
                          [&](uint i, uint j, uint k) {
        if (vox::isInsideSdf(boundarySdf.sample(pos(openvdb::Coord(i, j, k))))) {
            _markers(i, j, k) = kBoundary;
        } else if (vox::isInsideSdf(fluidSdf.sample(pos(openvdb::Coord(i, j, k))))) {
            _markers(i, j, k) = kFluid;
        } else {
            _markers(i, j, k) = kAir;
        }
    });
}
