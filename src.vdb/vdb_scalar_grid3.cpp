//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "../src.common/pch.h"

#include "vdb_fdm_utils.h"
#include "../src.common/parallel.h"
#include "vdb_scalar_grid3.h"
#include "vdb_samplers3.h"
#include "../src.common/serial.h"
#include <openvdb/math/Stencils.h>
#include <openvdb/math/Operators.h>
#include <algorithm>
#include <string>
#include <utility>  // just make cpplint happy..
#include <vector>

using namespace vdb;

ScalarGrid3::ScalarGrid3()
:_grid(openvdb::DoubleGrid::create(0.0)),
_linearSampler(LinearGridSampler<openvdb::DoubleGrid>(_grid,
                                                      vox::Size3(),
                                                      vox::Vector3D(1,1,1),
                                                      vox::Vector3D(0,0,0))) {}

ScalarGrid3::~ScalarGrid3() {}

void ScalarGrid3::clear() { resize(resolution(),
                                   gridSpacing(),
                                   origin(), 0.0); }

void ScalarGrid3::resize(size_t resolutionX,
                         size_t resolutionY,
                         size_t resolutionZ,
                         double gridSpacingX,
                         double gridSpacingY,
                         double gridSpacingZ,
                         double originX,
                         double originY,
                         double originZ,
                         double initialValue) {
    resize(vox::Size3(resolutionX, resolutionY, resolutionZ),
           vox::Vector3D(gridSpacingX, gridSpacingY, gridSpacingZ),
           vox::Vector3D(originX, originY, originZ),
           initialValue);
}

void ScalarGrid3::resize(const vox::Size3& resolution,
                         const vox::Vector3D& gridSpacing,
                         const vox::Vector3D& origin,
                         double initialValue) {
    setSizeParameters(resolution, gridSpacing, origin);
    _grid = openvdb::DoubleGrid::create(initialValue);
    resetSampler();
}

void ScalarGrid3::resize(double gridSpacingX,
                         double gridSpacingY,
                         double gridSpacingZ,
                         double originX,
                         double originY,
                         double originZ) {
    resize(vox::Vector3D(gridSpacingX, gridSpacingY, gridSpacingZ),
           vox::Vector3D(originX, originY, originZ));
}

void ScalarGrid3::resize(const vox::Vector3D& gridSpacing,
                         const vox::Vector3D& origin) {
    resize(resolution(), gridSpacing, origin);
}

double ScalarGrid3::operator()(const openvdb::Coord& coord) const {
    return _grid->tree().getValue(coord);
}

vox::Vector3D ScalarGrid3::gradientAtDataPoint(const openvdb::Coord& coord) const {
    return gradient3(_grid, dataSize(), gridSpacing(),
                     coord.x(), coord.y(), coord.z());
}

double ScalarGrid3::laplacianAtDataPoint(const openvdb::Coord& coord) const {
    return laplacian3(_grid, dataSize(), gridSpacing(),
                      coord.x(), coord.y(), coord.z());
}

double ScalarGrid3::sample(const vox::Vector3D& x) const { return _sampler(x); }

std::function<double(const vox::Vector3D&)> ScalarGrid3::sampler() const {
    return _sampler;
}

vox::Vector3D ScalarGrid3::gradient(const vox::Vector3D& x) const {
    std::array<openvdb::Coord, 8> indices;
    std::array<double, 8> weights;
    _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);
    
    vox::Vector3D result;
    
    for (int i = 0; i < 8; ++i) {
        result += weights[i] *
        gradientAtDataPoint(indices[i]);
    }
    
    return result;
}

double ScalarGrid3::laplacian(const vox::Vector3D& x) const {
    std::array<openvdb::Coord, 8> indices;
    std::array<double, 8> weights;
    _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);
    
    double result = 0.0;
    
    for (int i = 0; i < 8; ++i) {
        result += weights[i] * laplacianAtDataPoint(indices[i]);
    }
    
    return result;
}

openvdb::DoubleGrid::Accessor ScalarGrid3::dataAccessor(){
    return _grid->getAccessor();
}

openvdb::DoubleGrid::ConstAccessor ScalarGrid3::constDataAccessor() const{
    return _grid->getConstAccessor();
}

ScalarGrid3::DataPositionFunc ScalarGrid3::dataPosition() const {
    vox::Vector3D o = dataOrigin();
    return [this, o](const openvdb::Coord& coord) -> vox::Vector3D {
        return o + gridSpacing() * vox::Vector3D(coord.x(),
                                                 coord.y(),
                                                 coord.z());
    };
}

void ScalarGrid3::fill(double value, vox::ExecutionPolicy policy) {
    vox::parallelFor(
                     vox::kZeroSize, dataSize().x,
                     vox::kZeroSize, dataSize().y,
                     vox::kZeroSize, dataSize().z,
                     [this, value](uint i, uint j, uint k) {
        _grid->tree().setValueOnly(openvdb::Coord(i, j, k), value);
    }, policy);
}

void ScalarGrid3::fill(const std::function<double(const vox::Vector3D&)>& func,
                       vox::ExecutionPolicy policy) {
    DataPositionFunc pos = dataPosition();
    vox::parallelFor(vox::kZeroSize, dataSize().x,
                     vox::kZeroSize, dataSize().y,
                     vox::kZeroSize, dataSize().z,
                     [this, &func, &pos](uint i, uint j, uint k) {
        _grid->tree().setValueOnly(openvdb::Coord(i, j, k), func(pos(openvdb::Coord(i, j, k))));
    }, policy);
}

void ScalarGrid3::forEachDataPointIndex(
                                        const std::function<void(const openvdb::Coord&)>& func) const {
    for (openvdb::DoubleGrid::ValueOnIter iter = _grid->beginValueOn(); iter; ++iter) {
        func(iter.getCoord());
    }
}

void ScalarGrid3::parallelForEachDataPointIndex(
                                                const std::function<void(const openvdb::Coord&)>& func) const {
    using IterRange = openvdb::tree::IteratorRange<openvdb::DoubleGrid::ValueOnIter>;
    IterRange range(_grid->beginValueOn());
    tbb::parallel_for(range, [&func](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            openvdb::DoubleGrid::ValueOnIter iter = r.iterator();
            func(iter.getCoord());
        }
    });
}

void ScalarGrid3::swapScalarGrid(ScalarGrid3* other) {
    swapGrid(other);
    
    _grid.swap(other->_grid);
    std::swap(_linearSampler, other->_linearSampler);
    std::swap(_sampler, other->_sampler);
}

void ScalarGrid3::setScalarGrid(const ScalarGrid3& other) {
    setGrid(other);
    
    _grid = other._grid->deepCopy();
    resetSampler();
}

void ScalarGrid3::resetSampler() {
    _linearSampler = LinearGridSampler<openvdb::DoubleGrid>(_grid,
                                                            dataSize(),
                                                            gridSpacing(),
                                                            dataOrigin());
    _sampler = _linearSampler.functor();
}


void ScalarGrid3::getData(vox::ArrayAccessor3<double> data) const {
    JET_ASSERT(dataSize() == data.size());
    vox::serialFor(
                   vox::kZeroSize, dataSize().x,
                   vox::kZeroSize, dataSize().y,
                   vox::kZeroSize, dataSize().z,
                   [&](uint i, uint j, uint k) {
        data(i, j, k) = _grid->tree().getValue(openvdb::Coord(i, j, k) );
    });
}


void ScalarGrid3::setData(const vox::ConstArrayAccessor3<double> data) {
    JET_ASSERT(dataSize() == data.size());
    vox::serialFor(
                   vox::kZeroSize, dataSize().x,
                   vox::kZeroSize, dataSize().y,
                   vox::kZeroSize, dataSize().z,
                   [&](uint i, uint j, uint k) {
        _grid->tree().setValueOnly(openvdb::Coord(i, j, k), data(i, j, k) );
    });
}

ScalarGridBuilder3::ScalarGridBuilder3() {}

ScalarGridBuilder3::~ScalarGridBuilder3() {}
