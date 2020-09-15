//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"

#include "vdb_fmm_level_set_solver3.hpp"
#include "vdb_sph_points_to_implicit3.h"
#include "vdb_sph_system_data3.hpp"

using namespace vdb;

SphPointsToImplicit3::SphPointsToImplicit3(double kernelRadius,
                                           double cutOffDensity,
                                           bool isOutputSdf)
: _kernelRadius(kernelRadius),
_cutOffDensity(cutOffDensity),
_isOutputSdf(isOutputSdf) {}

void SphPointsToImplicit3::convert(const vox::ConstArrayAccessor1<vox::Vector3D>& points,
                                   ScalarGrid3* output) const {
    if (output == nullptr) {
        LOG(WARNING) << "Null scalar grid output pointer provided.";
        return;
    }
    
    const auto res = output->resolution();
    if (res.x * res.y * res.z == 0) {
        LOG(WARNING) << "Empty grid is provided.";
        return;
    }
    
    const auto bbox = output->boundingBox();
    if (bbox.isEmpty()) {
        LOG(WARNING) << "Empty domain is provided.";
        return;
    }
    
    SphSystemData3 sphParticles;
    sphParticles.addParticles(points);
    sphParticles.setKernelRadius(_kernelRadius);
    sphParticles.buildNeighborSearcher();
    sphParticles.updateDensities();
    
    vox::Array1<double> constData(sphParticles.numberOfParticles(), 1.0);
    auto temp = output->clone();
    temp->fill([&](const vox::Vector3D& x) {
        double d = sphParticles.interpolate(x, constData);
        return _cutOffDensity - d;
    }, vox::ExecutionPolicy::kSerial);
    
    if (_isOutputSdf) {
        FmmLevelSetSolver3 solver;
        solver.reinitialize(*temp, vox::kMaxD, output);
    } else {
        temp->swap(output);
    }
}
