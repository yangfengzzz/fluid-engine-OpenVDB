//
//  point_VDB_searcher3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/1/5.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"

#include "vdb_point_searcher3.hpp"

#include <algorithm>
#include <vector>

using namespace vdb;

PointVDBSearcher3::PointVDBSearcher3() {
}

PointVDBSearcher3::PointVDBSearcher3(
                                     const PointVDBSearcher3& other) {
    set(other);
}

//-----------------------------------------------------------------------------
void PointVDBSearcher3::build(const vox::ConstArrayAccessor1<vox::Vector3D>& points) {
    _points.resize(points.size());
    vox::parallelFor(
                     vox::kZeroSize,
                     points.size(),
                     [&](size_t i) {
        _points[i].x() = points[i].x;
        _points[i].y() = points[i].y;
        _points[i].z() = points[i].z;
    });
    
    openvdb::points::PointAttributeVector<openvdb::Vec3R> positionsWrapper(_points);
    
    int pointsPerVoxel = 8;
    float voxelSize = openvdb::points::computeVoxelSize(positionsWrapper, pointsPerVoxel);
    // Create a transform using this voxel-size.
    transform = openvdb::math::Transform::createLinearTransform(voxelSize);
    
    pointIndexGrid = openvdb::tools::createPointIndexGrid<openvdb::tools::PointIndexGrid>(positionsWrapper, *transform);
}

template<typename T>
struct Filter {
    typedef T ValueType;
    Filter(T const * const array, vox::PointNeighborSearcher3::ForEachNearbyPointFunc callback)
    : mValues(array), func(callback){}
    
    void operator()(const double distSqr, const size_t pointIndex) {
        func(pointIndex, vox::Vector3D(mValues[pointIndex].x(),
                                       mValues[pointIndex].y(),
                                       mValues[pointIndex].z()) );
    }
    
private:
    T const * const mValues;
    vox::PointNeighborSearcher3::ForEachNearbyPointFunc func;
}; // struct WeightedAverageAccumulator

void PointVDBSearcher3::forEachNearbyPoint(const vox::Vector3D& origin, double radius,
                                           const ForEachNearbyPointFunc& callback) {
    PointList pointList(_points);
    Filter<openvdb::Vec3d> filter(_points.data(), callback);
    openvdb::tools::PointIndexFilter<PointList> search(pointList, pointIndexGrid->tree(), *transform);
    search.searchAndApply(openvdb::Vec3d(origin.x, origin.y, origin.z),
                          radius, filter);
}

bool PointVDBSearcher3::hasNearbyPoint(const vox::Vector3D& origin, double radius) {
    PointList pointList(_points);
    
    openvdb::tools::PointIndexGrid::ConstAccessor acc = pointIndexGrid->getConstAccessor();
    
    iterator.worldSpaceSearchAndUpdate(openvdb::Vec3d(origin.x, origin.y, origin.z),
                                       radius, acc, pointList, *transform);
    
    if (iterator.size() != 0) {
        return true;
    }
    
    return false;
}

openvdb::tools::PointIndexGrid::Ptr PointVDBSearcher3::getIndexGrid(){
    return pointIndexGrid;
}

openvdb::points::PointDataGrid::Ptr PointVDBSearcher3::getDataGrid(){
    openvdb::points::PointAttributeVector<openvdb::Vec3R> positionsWrapper(_points);
    
    pointDataGrid =
    openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
    openvdb::points::PointDataGrid>(*pointIndexGrid, positionsWrapper, *transform);
    
    return pointDataGrid;
}

openvdb::math::Transform::Ptr PointVDBSearcher3::getTransform(){
    return transform;
}

openvdb::Vec3d PointVDBSearcher3::getPoint(size_t i){
    return _points[i];
}
//-----------------------------------------------------------------------------
vox::PointNeighborSearcher3Ptr PointVDBSearcher3::clone() const {
    return CLONE_W_CUSTOM_DELETER(PointVDBSearcher3);
}

PointVDBSearcher3&
PointVDBSearcher3::operator=(const PointVDBSearcher3& other) {
    set(other);
    return *this;
}

void PointVDBSearcher3::set(const PointVDBSearcher3& other) {
    _points = other._points;
    transform = other.transform;
    pointIndexGrid = other.pointIndexGrid;
}

//! Serializes the neighbor searcher into the buffer.
void PointVDBSearcher3::serialize(std::vector<uint8_t>* buffer) const{
    ;
}

//! Deserializes the neighbor searcher from the buffer.
void PointVDBSearcher3::deserialize(const std::vector<uint8_t>& buffer){
    ;
}

PointVDBSearcher3
PointVDBSearcher3::Builder::build() const {
    return PointVDBSearcher3();
}

PointVDBSearcher3Ptr
PointVDBSearcher3::Builder::makeShared() const {
    return std::shared_ptr<PointVDBSearcher3>(
                                              new PointVDBSearcher3(),
                                              [] (PointVDBSearcher3* obj) {
        delete obj;
    });
}

vox::PointNeighborSearcher3Ptr
PointVDBSearcher3::Builder::buildPointNeighborSearcher() const {
    return makeShared();
}
