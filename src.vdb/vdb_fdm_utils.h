//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_FDM_UTILS_H_
#define INCLUDE_VDB_FDM_UTILS_H_

#include <openvdb/openvdb.h>
#include "../src.common/vector3.h"
#include "../src.common/size3.h"
#include <iostream>

namespace vdb {
//! \brief Returns 3-D gradient vector from given 3-D scalar grid-like array
//!        \p data, \p gridSpacing, and array index (\p i, \p j, \p k).
vox::Vector3D gradient3(
                         const openvdb::DoubleGrid::Ptr data,
                         const vox::Size3 size,
                         const vox::Vector3D& gridSpacing,
                         uint i,
                         uint j,
                         uint k);

//! \brief Returns 3-D gradient vectors from given 3-D vector grid-like array
//!        \p data, \p gridSpacing, and array index (\p i, \p j, \p k).
std::array<vox::Vector3D, 3> gradient3(
                                        const openvdb::Vec3dGrid::Ptr data,
                                        const vox::Size3 size,
                                        const vox::Vector3D& gridSpacing,
                                        uint i,
                                        uint j,
                                        uint k);


//! \brief Returns Laplacian value from given 3-D scalar grid-like array
//!        \p data, \p gridSpacing, and array index (\p i, \p j, \p k).
double laplacian3(
                  const openvdb::DoubleGrid::Ptr data,
                  const vox::Size3 size,
                  const vox::Vector3D& gridSpacing,
                  uint i,
                  uint j,
                  uint k);

//! \brief Returns 3-D Laplacian vectors from given 3-D vector grid-like array
//!        \p data, \p gridSpacing, and array index (\p i, \p j, \p k).
vox::Vector3D laplacian3(
                          const openvdb::Vec3dGrid::Ptr data,
                          const vox::Size3 size,
                          const vox::Vector3D& gridSpacing,
                          uint i,
                          uint j,
                          uint k);


}  // namespace vox

#endif  // INCLUDE_VDB_FDM_UTILS_H_
