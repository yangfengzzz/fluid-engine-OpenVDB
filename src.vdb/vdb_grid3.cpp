//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"

#include "vdb_grid3.h"
#include "../src.common/parallel.h"
#include "../src.common/serial.h"
#include "../src.common/math_utils.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>  // just make cpplint happy..

using namespace vdb;

Grid3::Grid3() {}

Grid3::~Grid3() {}

const vox::Size3& Grid3::resolution() const { return _resolution; }

const vox::Vector3D& Grid3::origin() const { return _origin; }

const vox::Vector3D& Grid3::gridSpacing() const { return _gridSpacing; }

const vox::BoundingBox3D& Grid3::boundingBox() const { return _boundingBox; }

Grid3::DataPositionFunc Grid3::cellCenterPosition() const {
    vox::Vector3D h = _gridSpacing;
    vox::Vector3D o = _origin;
    return [h, o](const openvdb::Coord& coord) {
        return o + h * vox::Vector3D(coord.x() + 0.5,
                                     coord.y() + 0.5,
                                     coord.z() + 0.5);
    };
}

void Grid3::forEachCellIndex(
                             const std::function<void(uint, uint, uint)>& func) const {
    vox::serialFor(vox::kZeroSize, _resolution.x,
                   vox::kZeroSize, _resolution.y,
                   vox::kZeroSize, _resolution.z,
                   [&func](uint i, uint j, uint k) { func(i, j, k); });
}

void Grid3::parallelForEachCellIndex(
                                     const std::function<void(uint, uint, uint)>& func) const {
    vox::parallelFor(vox::kZeroSize, _resolution.x,
                     vox::kZeroSize, _resolution.y,
                     vox::kZeroSize, _resolution.z,
                     [&func](uint i, uint j, uint k) { func(i, j, k); });
}

bool Grid3::hasSameShape(const Grid3& other) const {
    return _resolution.x == other._resolution.x &&
    _resolution.y == other._resolution.y &&
    _resolution.z == other._resolution.z &&
    vox::similar(_gridSpacing.x, other._gridSpacing.x) &&
    vox::similar(_gridSpacing.y, other._gridSpacing.y) &&
    vox::similar(_gridSpacing.z, other._gridSpacing.z) &&
    vox::similar(_origin.x, other._origin.x) &&
    vox::similar(_origin.y, other._origin.y) &&
    vox::similar(_origin.z, other._origin.z);
}

void Grid3::setSizeParameters(const vox::Size3& resolution,
                              const vox::Vector3D& gridSpacing,
                              const vox::Vector3D& origin){
    _resolution = resolution;
    _origin = origin;
    _gridSpacing = gridSpacing;
    
    vox::Vector3D resolutionD = vox::Vector3D(static_cast<double>(resolution.x),
                                              static_cast<double>(resolution.y),
                                              static_cast<double>(resolution.z));
    
    _boundingBox = vox::BoundingBox3D(origin, origin + gridSpacing * resolutionD);
}

void Grid3::swapGrid(Grid3* other) {
    std::swap(_resolution, other->_resolution);
    std::swap(_gridSpacing, other->_gridSpacing);
    std::swap(_origin, other->_origin);
    std::swap(_boundingBox, other->_boundingBox);
}

void Grid3::setGrid(const Grid3& other) {
    _resolution = other._resolution;
    _gridSpacing = other._gridSpacing;
    _origin = other._origin;
    _boundingBox = other._boundingBox;
}
