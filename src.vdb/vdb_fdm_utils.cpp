//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"
#include "vdb_fdm_utils.h"

namespace vdb {
vox::Vector3D gradient3(
                         const openvdb::DoubleGrid::Ptr data,
                         const vox::Size3 ds,
                         const vox::Vector3D& gridSpacing,
                         uint i,
                         uint j,
                         uint k){
    JET_ASSERT(i < ds.x && j < ds.y && k < ds.z);
    
    double left = data->tree().getValue(openvdb::Coord((i > 0) ? i - 1 : i, j, k));
    double right = data->tree().getValue(openvdb::Coord((i + 1 < ds.x) ? i + 1 : i, j, k));
    double down = data->tree().getValue(openvdb::Coord(i, (j > 0) ? j - 1 : j, k));
    double up = data->tree().getValue(openvdb::Coord(i, (j + 1 < ds.y) ? j + 1 : j, k));
    double back = data->tree().getValue(openvdb::Coord(i, j, (k > 0) ? k - 1 : k));
    double front = data->tree().getValue(openvdb::Coord(i, j, (k + 1 < ds.z) ? k + 1 : k));
    
    return 0.5 * vox::Vector3D(right - left, up - down, front - back) / gridSpacing;
}

std::array<vox::Vector3D, 3> gradient3(
                                        const openvdb::Vec3dGrid::Ptr data,
                                        const vox::Size3 ds,
                                        const vox::Vector3D& gridSpacing,
                                        uint i,
                                        uint j,
                                        uint k){
    JET_ASSERT(i < ds.x && j < ds.y && k < ds.z);
    
    openvdb::Vec3d left = data->tree().getValue(openvdb::Coord((i > 0) ? i - 1 : i, j, k));
    openvdb::Vec3d right = data->tree().getValue(openvdb::Coord((i + 1 < ds.x) ? i + 1 : i, j, k));
    openvdb::Vec3d down = data->tree().getValue(openvdb::Coord(i, (j > 0) ? j - 1 : j, k));
    openvdb::Vec3d up = data->tree().getValue(openvdb::Coord(i, (j + 1 < ds.y) ? j + 1 : j, k));
    openvdb::Vec3d back = data->tree().getValue(openvdb::Coord(i, j, (k > 0) ? k - 1 : k));
    openvdb::Vec3d front = data->tree().getValue(openvdb::Coord(i, j, (k + 1 < ds.z) ? k + 1 : k));
    
    std::array<vox::Vector3D, 3> result;
    result[0] = 0.5 * vox::Vector3D(
                                     right.x() - left.x(), up.x() - down.x(), front.x() - back.x()) / gridSpacing;
    result[1] = 0.5 * vox::Vector3D(
                                     right.y() - left.y(), up.y() - down.y(), front.y() - back.y()) / gridSpacing;
    result[2] = 0.5 * vox::Vector3D(
                                     right.z() - left.z(), up.z() - down.z(), front.z() - back.z()) / gridSpacing;
    return result;
}


double laplacian3(
                  const openvdb::DoubleGrid::Ptr data,
                  const vox::Size3 ds,
                  const vox::Vector3D& gridSpacing,
                  uint i,
                  uint j,
                  uint k){
    const double center = data->tree().getValue(openvdb::Coord(i, j, k));

    JET_ASSERT(i < ds.x && j < ds.y && k < ds.z);
    
    double dleft = 0.0;
    double dright = 0.0;
    double ddown = 0.0;
    double dup = 0.0;
    double dback = 0.0;
    double dfront = 0.0;
    
    if (i > 0) {
        dleft = center - data->tree().getValue(openvdb::Coord(i - 1, j, k));
    }
    if (i + 1 < ds.x) {
        dright = data->tree().getValue(openvdb::Coord(i + 1, j, k)) - center;
    }

    if (j > 0) {
        ddown = center - data->tree().getValue(openvdb::Coord(i, j - 1, k));
    }
    if (j + 1 < ds.y) {
        dup = data->tree().getValue(openvdb::Coord(i, j + 1, k)) - center;
    }

    if (k > 0) {
        dback = center - data->tree().getValue(openvdb::Coord(i, j, k - 1));
    }
    if (k + 1 < ds.z) {
        dfront = data->tree().getValue(openvdb::Coord(i, j, k + 1)) - center;
    }
    
    return (dright - dleft) / vox::square(gridSpacing.x)
    + (dup - ddown) / vox::square(gridSpacing.y)
    + (dfront - dback) / vox::square(gridSpacing.z);
}

vox::Vector3D laplacian3(
                          const openvdb::Vec3dGrid::Ptr data,
                          const vox::Size3 ds,
                          const vox::Vector3D& gridSpacing,
                          uint i,
                          uint j,
                          uint k){
    const openvdb::Vec3d center = data->tree().getValue(openvdb::Coord(i, j, k));

    JET_ASSERT(i < ds.x && j < ds.y && k < ds.z);
    
    openvdb::Vec3d dleft = openvdb::Vec3d::zero();
    openvdb::Vec3d dright = openvdb::Vec3d::zero();
    openvdb::Vec3d ddown = openvdb::Vec3d::zero();
    openvdb::Vec3d dup = openvdb::Vec3d::zero();
    openvdb::Vec3d dback = openvdb::Vec3d::zero();
    openvdb::Vec3d dfront = openvdb::Vec3d::zero();
    
    if (i > 0) {
        dleft = center - data->tree().getValue(openvdb::Coord(i - 1, j, k));
    }
    if (i + 1 < ds.x) {
        dright = data->tree().getValue(openvdb::Coord(i + 1, j, k)) - center;
    }

    if (j > 0) {
        ddown = center - data->tree().getValue(openvdb::Coord(i, j - 1, k));
    }
    if (j + 1 < ds.y) {
        dup = data->tree().getValue(openvdb::Coord(i, j + 1, k)) - center;
    }

    if (k > 0) {
        dback = center - data->tree().getValue(openvdb::Coord(i, j, k - 1));
    }
    if (k + 1 < ds.z) {
        dfront = data->tree().getValue(openvdb::Coord(i, j, k + 1)) - center;
    }
    
    openvdb::Vec3d val =  (dright - dleft) / vox::square(gridSpacing.x)
    + (dup - ddown) / vox::square(gridSpacing.y)
    + (dfront - dback) / vox::square(gridSpacing.z);
    
    return vox::Vector3D(val.x(), val.y(), val.z());
}


}  // namespace vox
