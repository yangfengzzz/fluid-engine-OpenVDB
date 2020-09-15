// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef SRC_TESTS_UNIT_TESTS_UNIT_TESTS_UTILS_H_
#define SRC_TESTS_UNIT_TESTS_UNIT_TESTS_UTILS_H_

#include "../src.common/vector2.h"
#include <openvdb/openvdb.h>
#include "../external/gtest/include/gtest/gtest.h"

#define EXPECT_VECTOR2_EQ(expected, actual)     \
    EXPECT_DOUBLE_EQ((expected).x, (actual).x); \
    EXPECT_DOUBLE_EQ((expected).y, (actual).y);

#define EXPECT_VECTOR2_NEAR(expected, actual, eps) \
    EXPECT_NEAR((expected).x, (actual).x, eps);    \
    EXPECT_NEAR((expected).y, (actual).y, eps);

#define EXPECT_VECTOR3_EQ(expected, actual)     \
    EXPECT_DOUBLE_EQ((expected).x(), (actual).x()); \
    EXPECT_DOUBLE_EQ((expected).y(), (actual).y()); \
    EXPECT_DOUBLE_EQ((expected).z(), (actual).z());

#define EXPECT_VECTOR3_NEAR(expected, actual, eps) \
    EXPECT_NEAR((expected).x(), (actual).x(), eps);    \
    EXPECT_NEAR((expected).y(), (actual).y(), eps);    \
    EXPECT_NEAR((expected).z(), (actual).z(), eps);

#define EXPECT_VECTOR4_EQ(expected, actual)     \
    EXPECT_DOUBLE_EQ((expected).x(), (actual).x()); \
    EXPECT_DOUBLE_EQ((expected).y(), (actual).y()); \
    EXPECT_DOUBLE_EQ((expected).z(), (actual).z()); \
    EXPECT_DOUBLE_EQ((expected).w(), (actual).w());

#define EXPECT_VECTOR4_NEAR(expected, actual, eps) \
    EXPECT_NEAR((expected).x(), (actual).x(), eps);    \
    EXPECT_NEAR((expected).y(), (actual).y(), eps);    \
    EXPECT_NEAR((expected).z(), (actual).z(), eps);    \
    EXPECT_NEAR((expected).w(), (actual).w(), eps);

#define EXPECT_BOUNDING_BOX2_EQ(expected, actual)                    \
    EXPECT_VECTOR2_EQ((expected).lowerCorner, (actual).lowerCorner); \
    EXPECT_VECTOR2_EQ((expected).upperCorner, (actual).upperCorner);

#define EXPECT_BOUNDING_BOX2_NEAR(expected, actual, eps)                    \
    EXPECT_VECTOR2_NEAR((expected).lowerCorner, (actual).lowerCorner, eps); \
    EXPECT_VECTOR2_NEAR((expected).upperCorner, (actual).upperCorner, eps);

#define EXPECT_BOUNDING_BOX3_EQ(expected, actual)                    \
    EXPECT_VECTOR3_EQ((expected).min(), (actual).min()); \
    EXPECT_VECTOR3_EQ((expected).max(), (actual).max());

#define EXPECT_BOUNDING_BOX3_NEAR(expected, actual, eps)                    \
    EXPECT_VECTOR3_NEAR((expected).min(), (actual).min(), eps); \
    EXPECT_VECTOR3_NEAR((expected).max(), (actual).max(), eps);

namespace vdb {

const vox::Vector2D* getSamplePoints2();

size_t getNumberOfSamplePoints2();

const openvdb::Vec3d* getSamplePoints3();

size_t getNumberOfSamplePoints3();

const vox::Vector2D* getSampleDirs2();

size_t getNumberOfSampleDirs2();

const openvdb::Vec3d* getSampleDirs3();

size_t getNumberOfSampleDirs3();

const char* getCubeTriMesh3x3x3Obj();

const char* getSphereTriMesh5x5Obj();

openvdb::Vec3d curlAtDataPoint(openvdb::Vec3DGrid::Ptr grid,
                               size_t i, size_t j, size_t k);

double divergenceAtDataPoint(openvdb::Vec3DGrid::Ptr grid,
                             size_t i, size_t j, size_t k);

}  // namespace vox

#endif  // SRC_TESTS_UNIT_TESTS_UNIT_TESTS_UTILS_H_
