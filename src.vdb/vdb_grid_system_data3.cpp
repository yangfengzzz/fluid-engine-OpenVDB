//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "../src.common/pch.h"

#include "vdb_grid_system_data3.h"

#include <algorithm>
#include <vector>

using namespace vdb;

GridSystemData3::GridSystemData3()
: GridSystemData3({0, 0, 0}, {1, 1, 1}, {0, 0, 0}) {
}

GridSystemData3::GridSystemData3(const vox::Size3& resolution,
                                 const vox::Vector3D& gridSpacing,
                                 const vox::Vector3D& origin) {
    _velocity = std::make_shared<FaceCenteredGrid3>();
    _advectableVectorDataList.push_back(_velocity);
    _velocityIdx = 0;
    resize(resolution, gridSpacing, origin);
}

GridSystemData3::GridSystemData3(const GridSystemData3& other) {
    resize(other._resolution, other._gridSpacing, other._origin);
    
    for (auto& data : other._scalarDataList) {
        _scalarDataList.push_back(data->clone());
    }
    for (auto& data : other._vectorDataList) {
        _vectorDataList.push_back(data->clone());
    }
    for (auto& data : other._advectableScalarDataList) {
        _advectableScalarDataList.push_back(data->clone());
    }
    for (auto& data : other._advectableVectorDataList) {
        _advectableVectorDataList.push_back(data->clone());
    }
    
    JET_ASSERT(_advectableVectorDataList.size() > 0);
    
    _velocity = std::dynamic_pointer_cast<FaceCenteredGrid3>(
                                                             _advectableVectorDataList[0]);
    
    JET_ASSERT(_velocity != nullptr);
    
    _velocityIdx = 0;
}

GridSystemData3::~GridSystemData3() {
}

void GridSystemData3::resize(
                             const vox::Size3& resolution,
                             const vox::Vector3D& gridSpacing,
                             const vox::Vector3D& origin) {
    _resolution = resolution;
    _gridSpacing = gridSpacing;
    _origin = origin;
    
    for (auto& data : _scalarDataList) {
        data->resize(resolution, gridSpacing, origin);
    }
    for (auto& data : _vectorDataList) {
        data->resize(resolution, gridSpacing, origin);
    }
    for (auto& data : _advectableScalarDataList) {
        data->resize(resolution, gridSpacing, origin);
    }
    for (auto& data : _advectableVectorDataList) {
        data->resize(resolution, gridSpacing, origin);
    }
}

vox::Size3 GridSystemData3::resolution() const {
    return _resolution;
}

vox::Vector3D GridSystemData3::gridSpacing() const {
    return _gridSpacing;
}

vox::Vector3D GridSystemData3::origin() const {
    return _origin;
}

vox::BoundingBox3D GridSystemData3::boundingBox() const {
    return _velocity->boundingBox();
}

size_t GridSystemData3::addScalarData(
                                      const ScalarGridBuilder3Ptr& builder,
                                      double initialVal) {
    size_t attrIdx = _scalarDataList.size();
    _scalarDataList.push_back(
                              builder->build(resolution(),
                                             gridSpacing(),
                                             origin(),
                                             initialVal));
    return attrIdx;
}

size_t GridSystemData3::addVectorData(
                                      const VectorGridBuilder3Ptr& builder,
                                      const vox::Vector3D& initialVal) {
    size_t attrIdx = _vectorDataList.size();
    _vectorDataList.push_back(
                              builder->build(resolution(),
                                             gridSpacing(),
                                             origin(),
                                             initialVal));
    return attrIdx;
}

size_t GridSystemData3::addAdvectableScalarData(
                                                const ScalarGridBuilder3Ptr& builder,
                                                double initialVal) {
    size_t attrIdx = _advectableScalarDataList.size();
    _advectableScalarDataList.push_back(
                                        builder->build(resolution(),
                                                       gridSpacing(),
                                                       origin(),
                                                       initialVal));
    return attrIdx;
}

size_t GridSystemData3::addAdvectableVectorData(
                                                const VectorGridBuilder3Ptr& builder,
                                                const vox::Vector3D& initialVal) {
    size_t attrIdx = _advectableVectorDataList.size();
    _advectableVectorDataList.push_back(
                                        builder->build(resolution(),
                                                       gridSpacing(),
                                                       origin(),
                                                       initialVal));
    return attrIdx;
}

const FaceCenteredGrid3Ptr& GridSystemData3::velocity() const {
    return _velocity;
}

size_t GridSystemData3::velocityIndex() const {
    return _velocityIdx;
}

const ScalarGrid3Ptr& GridSystemData3::scalarDataAt(size_t idx) const {
    return _scalarDataList[idx];
}

const VectorGrid3Ptr& GridSystemData3::vectorDataAt(size_t idx) const {
    return _vectorDataList[idx];
}

const ScalarGrid3Ptr&
GridSystemData3::advectableScalarDataAt(size_t idx) const {
    return _advectableScalarDataList[idx];
}

const VectorGrid3Ptr&
GridSystemData3::advectableVectorDataAt(size_t idx) const {
    return _advectableVectorDataList[idx];
}

size_t GridSystemData3::numberOfScalarData() const {
    return _scalarDataList.size();
}

size_t GridSystemData3::numberOfVectorData() const {
    return _vectorDataList.size();
}

size_t GridSystemData3::numberOfAdvectableScalarData() const {
    return _advectableScalarDataList.size();
}

size_t GridSystemData3::numberOfAdvectableVectorData() const {
    return _advectableVectorDataList.size();
}
