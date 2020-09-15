//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"

#include "vdb_fmm_level_set_solver3.hpp"
#include "vdb_particle_system_data3.hpp"
#include "vdb_zhu_bridson_points_to_implicit3.h"

using namespace vdb;

inline double k(double s) { return std::max(0.0, vox::cubic(1.0 - s * s)); }

ZhuBridsonPointsToImplicit3::ZhuBridsonPointsToImplicit3(double kernelRadius,
                                                         double cutOffThreshold,
                                                         bool isOutputSdf)
: _kernelRadius(kernelRadius),
_cutOffThreshold(cutOffThreshold),
_isOutputSdf(isOutputSdf) {}

void ZhuBridsonPointsToImplicit3::convert(const vox::ConstArrayAccessor1<vox::Vector3D>& points,
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
    
    ParticleSystemData3 particles;
    particles.addParticles(points);
    particles.buildNeighborSearcher(_kernelRadius);
    
    const auto neighborSearcher = particles.neighborSearcher();
    const double isoContValue = _cutOffThreshold * _kernelRadius;
    
    auto temp = output->clone();
    openvdb::CoordBBox box = output->getGrid()->evalActiveVoxelBoundingBox();
    double L= output->getGrid()->transform().indexToWorld(box).extents().length();
    temp->fill([&](const vox::Vector3D& x) -> double {
        vox::Vector3D xAvg;
        double wSum = 0.0;
        const auto func = [&](size_t, const vox::Vector3D& xi) {
            const double wi = k((x - xi).length() / _kernelRadius);
            wSum += wi;
            xAvg += wi * xi;
        };
        neighborSearcher->forEachNearbyPoint(x, _kernelRadius, func);
        
        if (wSum > 0.0) {
            xAvg /= wSum;
            return (x - xAvg).length() - isoContValue;
        } else {
            return L;
        }
    }, vox::ExecutionPolicy::kSerial);
    
    if (_isOutputSdf) {
        FmmLevelSetSolver3 solver;
        solver.reinitialize(*temp, vox::kMaxD, output);
    } else {
        temp->swap(output);
    }
}
