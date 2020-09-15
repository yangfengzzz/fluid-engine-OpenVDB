//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"
#include "vdb_collocated_vector_grid3.h"
#include "../src.common/parallel.h"
#include "../src.common/serial.h"
#include "vdb_samplers3.h"
#include <openvdb/math/Stencils.h>
#include <openvdb/math/Operators.h>
#include <algorithm>
#include <utility>  // just make cpplint happy..
#include <vector>

using namespace vdb;

CollocatedVectorGrid3::CollocatedVectorGrid3() :
_linearSampler(_grid, vox::Size3(),
               vox::Vector3D(1, 1, 1),
               vox::Vector3D()) {}

CollocatedVectorGrid3::~CollocatedVectorGrid3() {
}

openvdb::Vec3d CollocatedVectorGrid3::operator()(const openvdb::Coord& coord) {
    return _grid->tree().getValue(coord);
}

double CollocatedVectorGrid3::divergenceAtDataPoint(const openvdb::Coord& coord) const {
    const vox::Size3 ds = dataSize();
    const vox::Vector3D& gs = gridSpacing();
    int i = coord.x();
    int j = coord.y();
    int k = coord.z();
    
    JET_ASSERT(coord.x() < ds.x && coord.y() < ds.y && coord.z() < ds.z);
    
    double left = _grid->tree().getValue(openvdb::Coord((i > 0) ? i - 1 : i, j, k)).x();
    double right = _grid->tree().getValue(openvdb::Coord((i + 1 < ds.x) ? i + 1 : i, j, k)).x();
    double down = _grid->tree().getValue(openvdb::Coord(i, (j > 0) ? j - 1 : j, k)).y();
    double up = _grid->tree().getValue(openvdb::Coord(i, (j + 1 < ds.y) ? j + 1 : j, k)).y();
    double back = _grid->tree().getValue(openvdb::Coord(i, j, (k > 0) ? k - 1 : k)).z();
    double front = _grid->tree().getValue(openvdb::Coord(i, j, (k + 1 < ds.z) ? k + 1 : k)).z();
    
    return 0.5 * (right - left) / gs.x
    + 0.5 * (up - down) / gs.y
    + 0.5 * (front - back) / gs.z;
}

vox::Vector3D CollocatedVectorGrid3::curlAtDataPoint(const openvdb::Coord& coord) const {
    const vox::Size3 ds = dataSize();
    const vox::Vector3D& gs = gridSpacing();
    int i = coord.x();
    int j = coord.y();
    int k = coord.z();
    
    JET_ASSERT(i < ds.x && j < ds.y && k < ds.z);
    
    openvdb::Vec3d left = _grid->tree().getValue(openvdb::Coord((i > 0) ? i - 1 : i, j, k));
    openvdb::Vec3d right = _grid->tree().getValue(openvdb::Coord((i + 1 < ds.x) ? i + 1 : i, j, k));
    openvdb::Vec3d down = _grid->tree().getValue(openvdb::Coord(i, (j > 0) ? j - 1 : j, k));
    openvdb::Vec3d up = _grid->tree().getValue(openvdb::Coord(i, (j + 1 < ds.y) ? j + 1 : j, k));
    openvdb::Vec3d back = _grid->tree().getValue(openvdb::Coord(i, j, (k > 0) ? k - 1 : k));
    openvdb::Vec3d front = _grid->tree().getValue(openvdb::Coord(i, j, (k + 1 < ds.z) ? k + 1 : k));
    
    double Fx_ym = down.x();
    double Fx_yp = up.x();
    double Fx_zm = back.x();
    double Fx_zp = front.x();
    
    double Fy_xm = left.y();
    double Fy_xp = right.y();
    double Fy_zm = back.y();
    double Fy_zp = front.y();
    
    double Fz_xm = left.z();
    double Fz_xp = right.z();
    double Fz_ym = down.z();
    double Fz_yp = up.z();
    
    return vox::Vector3D(
                          0.5 * (Fz_yp - Fz_ym) / gs.y - 0.5 * (Fy_zp - Fy_zm) / gs.z,
                          0.5 * (Fx_zp - Fx_zm) / gs.z - 0.5 * (Fz_xp - Fz_xm) / gs.x,
                          0.5 * (Fy_xp - Fy_xm) / gs.x - 0.5 * (Fx_yp - Fx_ym) / gs.y);
}

vox::Vector3D CollocatedVectorGrid3::sample(const vox::Vector3D& x) const {
    return _sampler(x);
}

double CollocatedVectorGrid3::divergence(const vox::Vector3D& x) const {
    std::array<openvdb::Coord, 8> indices;
    std::array<double, 8> weights;
    _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);
    
    double result = 0.0;
    
    for (int i = 0; i < 8; ++i) {
        result += weights[i] * divergenceAtDataPoint(indices[i]);
    }
    
    return result;
}

vox::Vector3D CollocatedVectorGrid3::curl(const vox::Vector3D& x) const {
    std::array<openvdb::Coord, 8> indices;
    std::array<double, 8> weights;
    _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);
    
    vox::Vector3D result;
    
    for (int i = 0; i < 8; ++i) {
        result += weights[i] * curlAtDataPoint(indices[i]);
    }
    
    return result;
}

std::function<vox::Vector3D(const vox::Vector3D&)>
CollocatedVectorGrid3::sampler() const {
    return _sampler;
}

openvdb::Vec3dGrid::Accessor CollocatedVectorGrid3::dataAccessor(){
    return _grid->getAccessor();
}

openvdb::Vec3dGrid::ConstAccessor CollocatedVectorGrid3::constDataAccessor() const{
    return _grid->getConstAccessor();
}

VectorGrid3::DataPositionFunc CollocatedVectorGrid3::dataPosition() const {
    vox::Vector3D dataOrigin_ = dataOrigin();
    return [this, dataOrigin_](const openvdb::Coord& coord) -> vox::Vector3D {
        return dataOrigin_ + gridSpacing() * vox::Vector3D(coord.x(),
                                                            coord.y(),
                                                            coord.z());
    };
}

void CollocatedVectorGrid3::forEachDataPointIndex(
                                                  const std::function<void(const openvdb::Coord& coord)>& func) const {
    for (openvdb::Vec3dGrid::ValueOnIter iter = _grid->beginValueOn(); iter; ++iter) {
        func(iter.getCoord());
    }
}

void CollocatedVectorGrid3::parallelForEachDataPointIndex(
                                                          const std::function<void(const openvdb::Coord& coord)>& func) const {
    using IterRange = openvdb::tree::IteratorRange<openvdb::Vec3dGrid::ValueOnIter>;
    IterRange range(_grid->beginValueOn());
    tbb::parallel_for(range, [&func](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            openvdb::Vec3dGrid::ValueOnIter iter = r.iterator();
            func(iter.getCoord());
        }
    });
}

void CollocatedVectorGrid3::fill(const vox::Vector3D& value,
                                 vox::ExecutionPolicy policy) {
    vox::parallelFor(
                     vox::kZeroSize, dataSize().x,
                     vox::kZeroSize, dataSize().y,
                     vox::kZeroSize, dataSize().z,
                     [this, value](uint i, uint j, uint k) {
        _grid->tree().setValueOnly(openvdb::Coord(i, j, k),
                                   openvdb::Vec3d(value.x,
                                                  value.y,
                                                  value.z) );
    }, policy);
}

void CollocatedVectorGrid3::fill(
                                 const std::function<vox::Vector3D(const vox::Vector3D&)>& func,
                                 vox::ExecutionPolicy policy) {
    DataPositionFunc pos = dataPosition();
    vox::parallelFor(vox::kZeroSize, dataSize().x,
                     vox::kZeroSize, dataSize().y,
                     vox::kZeroSize, dataSize().z,
                [this, &func, &pos](uint i, uint j, uint k) {
        vox::Vector3D value = func(pos(openvdb::Coord(i, j, k)));
        _grid->tree().setValueOnly(openvdb::Coord(i, j, k),
                                   openvdb::Vec3d(value.x,
                                                  value.y,
                                                  value.z) );
    }, policy);
}

void CollocatedVectorGrid3::swapCollocatedVectorGrid(
                                                     CollocatedVectorGrid3* other) {
    swapGrid(other);
    
    _grid.swap(other->_grid);
    std::swap(_linearSampler, other->_linearSampler);
    std::swap(_sampler, other->_sampler);
}

void CollocatedVectorGrid3::setCollocatedVectorGrid(
                                                    const CollocatedVectorGrid3& other) {
    setGrid(other);
    
    _grid = other._grid->deepCopy();
    resetSampler();
}

void CollocatedVectorGrid3::onResize(
                                     const vox::Size3& resolution,
                                     const vox::Vector3D& gridSpacing,
                                     const vox::Vector3D& origin,
                                     const vox::Vector3D& initialValue) {
    UNUSED_VARIABLE(resolution);
    UNUSED_VARIABLE(gridSpacing);
    UNUSED_VARIABLE(origin);
    _grid = openvdb::Vec3dGrid::create(openvdb::Vec3d(initialValue.x,
                                                      initialValue.y,
                                                      initialValue.z ));
    resetSampler();
}

void CollocatedVectorGrid3::resetSampler() {
    _linearSampler = LinearGridSampler<openvdb::Vec3dGrid>(_grid,
                                                           dataSize(),
                                                           gridSpacing(),
                                                           dataOrigin());
    _sampler = [&](const vox::Vector3D& pt)->vox::Vector3D{
        openvdb::Vec3d val = _linearSampler.functor()(pt);
        return vox::Vector3D(val.x(), val.y(), val.z());
    };
}
