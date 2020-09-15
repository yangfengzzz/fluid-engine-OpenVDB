//
//  vdb_cubic_semi_lagrangian3.cpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/21.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.common/pch.h"
#include "vdb_samplers3.h"
#include "vdb_cubic_semi_lagrangian3.hpp"

using namespace vdb;

CubicSemiLagrangian3::CubicSemiLagrangian3() {
}

std::function<double(const vox::Vector3D&)>
CubicSemiLagrangian3::getScalarSamplerFunc(const ScalarGrid3& source) const {
    auto sourceSampler = CubicGridSampler<openvdb::DoubleGrid>(source.getGrid(),
                                                               source.dataSize(),
                                                               source.gridSpacing(),
                                                               source.dataOrigin());
    return sourceSampler.functor();
}

std::function<vox::Vector3D(const vox::Vector3D&)>
CubicSemiLagrangian3::getVectorSamplerFunc(
                                           const CollocatedVectorGrid3& source) const {
    auto sourceSampler = CubicGridSampler<openvdb::Vec3dGrid>(
                                                              source.getGrid(),
                                                              source.dataSize(),
                                                              source.gridSpacing(),
                                                              source.dataOrigin());
    return [&](const vox::Vector3D& pt)->vox::Vector3D{
        openvdb::Vec3d val = sourceSampler(pt);
        return vox::Vector3D(val.x(), val.y(), val.z());
    };
}
