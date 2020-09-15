//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"

#include "vdb_fmm_level_set_solver3.hpp"
#include "vdb_particle_system_data3.hpp"
#include "vdb_spherical_points_to_implicit3.h"

using namespace vdb;

SphericalPointsToImplicit3::SphericalPointsToImplicit3(double radius,
                                                       bool isOutputSdf)
: _radius(radius), _isOutputSdf(isOutputSdf) {}

void SphericalPointsToImplicit3::convert(
                                         const vox::ConstArrayAccessor1<vox::Vector3D>& points,
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
    particles.buildNeighborSearcher(2.0 * _radius);
    
    const auto neighborSearcher = particles.neighborSearcher();
    
    auto temp = output->clone();
    temp->fill([&](const vox::Vector3D& x) {
        double minDist = 2.0 * _radius;
        neighborSearcher->forEachNearbyPoint(
                                             x,
                                             2.0 * _radius,
                                             [&](size_t, const vox::Vector3D& xj) {
            minDist = std::min(minDist, (x - xj).length());
        });
        
        return minDist - _radius;
    }, vox::ExecutionPolicy::kSerial);
    
    if (_isOutputSdf) {
        FmmLevelSetSolver3 solver;
        solver.reinitialize(*temp, vox::kMaxD, output);
    } else {
        temp->swap(output);
    }
}
