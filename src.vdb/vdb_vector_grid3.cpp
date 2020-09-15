//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "../src.common/pch.h"

#include "vdb_vector_grid3.h"

#include <algorithm>
#include <string>
#include <vector>

using namespace vdb;

VectorGrid3::VectorGrid3() {}

VectorGrid3::~VectorGrid3() {}

void VectorGrid3::clear() {
    resize(resolution(), gridSpacing(), origin(), vox::Vector3D());
}

void VectorGrid3::resize(size_t resolutionX,
                         size_t resolutionY,
                         size_t resolutionZ,
                         double gridSpacingX,
                         double gridSpacingY,
                         double gridSpacingZ,
                         double originX,
                         double originY,
                         double originZ,
                         double initialValueX,
                         double initialValueY,
                         double initialValueZ) {
    resize(vox::Size3(resolutionX, resolutionY, resolutionZ),
           vox::Vector3D(gridSpacingX, gridSpacingY, gridSpacingZ),
           vox::Vector3D(originX, originY, originZ),
           vox::Vector3D(initialValueX, initialValueY, initialValueZ));
}

void VectorGrid3::resize(const vox::Size3& resolution,
                         const vox::Vector3D& gridSpacing,
                         const vox::Vector3D& origin,
                         const vox::Vector3D& initialValue) {
    setSizeParameters(resolution, gridSpacing, origin);
    
    onResize(resolution, gridSpacing, origin, initialValue);
}

void VectorGrid3::resize(double gridSpacingX,
                         double gridSpacingY,
                         double gridSpacingZ,
                         double originX,
                         double originY,
                         double originZ) {
    resize(vox::Vector3D(gridSpacingX, gridSpacingY, gridSpacingZ),
           vox::Vector3D(originX, originY, originZ));
}

void VectorGrid3::resize(const vox::Vector3D& gridSpacing,
                         const vox::Vector3D& origin) {
    resize(resolution(), gridSpacing, origin);
}

VectorGridBuilder3::VectorGridBuilder3() {}

VectorGridBuilder3::~VectorGridBuilder3() {}
