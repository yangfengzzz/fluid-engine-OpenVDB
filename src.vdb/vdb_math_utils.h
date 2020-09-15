//
//  math_utils.h
//
//  Created by Feng Yang on 2020/5/15.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef math_utils_h
#define math_utils_h

#include <openvdb/openvdb.h>

namespace vox {

inline openvdb::Vec3d monotonicCatmullRom(const openvdb::Vec3d& v0, const openvdb::Vec3d& v1,
                                          const openvdb::Vec3d& v2, const openvdb::Vec3d& v3,
                                          double f) {
    static const double two = 2.0;
    static const double three = 3.0;

    openvdb::Vec3d d1 = (v2 - v0) / two;
    openvdb::Vec3d d2 = (v3 - v1) / two;
    openvdb::Vec3d D1 = v2 - v1;

    if (std::fabs(D1.x()) < std::numeric_limits<float>::epsilon() ||
        sign(D1.x()) != sign(d1.x()) || sign(D1.x()) != sign(d2.x())) {
        d1.x() = d2.x() = 0;
    }

    if (std::fabs(D1.y()) < std::numeric_limits<float>::epsilon() ||
        sign(D1.y()) != sign(d1.y()) || sign(D1.y()) != sign(d2.y())) {
        d1.y() = d2.y() = 0;
    }

    if (std::fabs(D1.z()) < std::numeric_limits<float>::epsilon() ||
        sign(D1.z()) != sign(d1.z()) || sign(D1.z()) != sign(d2.z())) {
        d1.z() = d2.z() = 0;
    }

    openvdb::Vec3d a3 = d1 + d2 - two * D1;
    openvdb::Vec3d a2 = three * D1 - two * d1 - d2;
    openvdb::Vec3d a1 = d1;
    openvdb::Vec3d a0 = v1;

    return a3 * cubic(f) + a2 * square(f) + a1 * f + a0;
}


}

#endif /* math_utils_h */
