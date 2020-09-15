//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"

#include "../src.common/parallel.h"
#include "vdb_vertex_centered_vector_grid3.h"

#include <utility>  // just make cpplint happy..

using namespace vdb;

VertexCenteredVectorGrid3::VertexCenteredVectorGrid3() :CollocatedVectorGrid3() {}

VertexCenteredVectorGrid3::VertexCenteredVectorGrid3(size_t resolutionX,
                                                     size_t resolutionY,
                                                     size_t resolutionZ,
                                                     double gridSpacingX,
                                                     double gridSpacingY,
                                                     double gridSpacingZ,
                                                     double originX,
                                                     double originY,
                                                     double originZ,
                                                     double initialValueU,
                                                     double initialValueV,
                                                     double initialValueW) {
    resize(resolutionX, resolutionY, resolutionZ, gridSpacingX, gridSpacingY,
           gridSpacingZ, originX, originY, originZ, initialValueU,
           initialValueV, initialValueW);
}

VertexCenteredVectorGrid3::VertexCenteredVectorGrid3(const vox::Size3& resolution,
                                                     const vox::Vector3D& gridSpacing,
                                                     const vox::Vector3D& origin,
                                                     const vox::Vector3D& initialValue) {
    resize(resolution, gridSpacing, origin, initialValue);
}

vox::Size3 VertexCenteredVectorGrid3::dataSize() const {
    if (resolution() != vox::Size3(0, 0, 0)) {
        return resolution() + vox::Size3(1, 1, 1);
    } else {
        return vox::Size3(0, 0, 0);
    }
}

vox::Vector3D VertexCenteredVectorGrid3::dataOrigin() const { return origin(); }

void VertexCenteredVectorGrid3::swap(Grid3* other) {
    VertexCenteredVectorGrid3* sameType =
    dynamic_cast<VertexCenteredVectorGrid3*>(other);
    if (sameType != nullptr) {
        swapCollocatedVectorGrid(sameType);
    }
}

std::shared_ptr<VectorGrid3> VertexCenteredVectorGrid3::clone() const {
    return CLONE_W_CUSTOM_DELETER(VertexCenteredVectorGrid3);
}

void VertexCenteredVectorGrid3::set(const VertexCenteredVectorGrid3& other) {
    setCollocatedVectorGrid(other);
}

VertexCenteredVectorGrid3& VertexCenteredVectorGrid3::operator=(
                                                                const VertexCenteredVectorGrid3& other) {
    set(other);
    return *this;
}

//------------------------------------------------------------------
VertexCenteredVectorGrid3::Builder VertexCenteredVectorGrid3::builder() {
    return Builder();
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withResolution(const vox::Size3& resolution) {
    _resolution = resolution;
    return *this;
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withResolution(size_t resolutionX,
                                                   size_t resolutionY,
                                                   size_t resolutionZ) {
    _resolution.x = resolutionX;
    _resolution.y = resolutionY;
    _resolution.z = resolutionZ;
    return *this;
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withGridSpacing(
                                                    const vox::Vector3D& gridSpacing) {
    _gridSpacing = gridSpacing;
    return *this;
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withGridSpacing(double gridSpacingX,
                                                    double gridSpacingY,
                                                    double gridSpacingZ) {
    _gridSpacing.x = gridSpacingX;
    _gridSpacing.y = gridSpacingY;
    _gridSpacing.z = gridSpacingZ;
    return *this;
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withOrigin(const vox::Vector3D& gridOrigin) {
    _gridOrigin = gridOrigin;
    return *this;
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withOrigin(double gridOriginX,
                                               double gridOriginY,
                                               double gridOriginZ) {
    _gridOrigin.x = gridOriginX;
    _gridOrigin.y = gridOriginY;
    _gridOrigin.z = gridOriginZ;
    return *this;
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withInitialValue(
                                                     const vox::Vector3D& initialVal) {
    _initialVal = initialVal;
    return *this;
}

VertexCenteredVectorGrid3::Builder&
VertexCenteredVectorGrid3::Builder::withInitialValue(double initialValX,
                                                     double initialValY,
                                                     double initialValZ) {
    _initialVal.x = initialValX;
    _initialVal.y = initialValY;
    _initialVal.z = initialValZ;
    return *this;
}

VertexCenteredVectorGrid3 VertexCenteredVectorGrid3::Builder::build() const {
    return VertexCenteredVectorGrid3(_resolution,
                                     _gridSpacing,
                                     _gridOrigin,
                                     _initialVal);
}

VertexCenteredVectorGrid3Ptr VertexCenteredVectorGrid3::Builder::makeShared()
const {
    return std::shared_ptr<VertexCenteredVectorGrid3>(
                                                      new VertexCenteredVectorGrid3(
                                                                                    _resolution,
                                                                                    _gridSpacing,
                                                                                    _gridOrigin,
                                                                                    _initialVal),
                                                      [](VertexCenteredVectorGrid3* obj) { delete obj; });
}

VectorGrid3Ptr VertexCenteredVectorGrid3::Builder::build(const vox::Size3& resolution,
                                                         const vox::Vector3D& gridSpacing,
                                                         const vox::Vector3D& gridOrigin,
                                                         const vox::Vector3D& initialVal) const {
    return std::shared_ptr<VertexCenteredVectorGrid3>(
                                                      new VertexCenteredVectorGrid3(
                                                                                    resolution,
                                                                                    gridSpacing,
                                                                                    gridOrigin,
                                                                                    initialVal),
                                                      [](VertexCenteredVectorGrid3* obj) { delete obj; });
}
