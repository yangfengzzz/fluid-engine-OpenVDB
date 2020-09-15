//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"

//#include "../src.common/math_utils.h"
#include "vdb_anisotropic_points_to_implicit3.h"
#include "vdb_fmm_level_set_solver3.hpp"
#include "../src.common/point_kdtree_searcher3.h"
#include "../src.common/sph_kernels3.h"
#include "vdb_sph_system_data3.hpp"
#include "../src.common/svd.h"

using namespace vdb;

inline double p(double distance) {
    const double distanceSquared = distance * distance;
    
    if (distanceSquared >= 1.0) {
        return 0.0;
    } else {
        const double x = 1.0 - distanceSquared;
        return x * x * x;
    }
}

inline double wij(double distance, double r) {
    if (distance < r) {
        return 1.0 - vox::cubic(distance / r);
    } else {
        return 0.0;
    }
}

inline vox::Matrix3x3D vvt(const vox::Vector3D& v) {
    return vox::Matrix3x3D(v.x * v.x, v.x * v.y, v.x * v.z,
                           v.y * v.x, v.y * v.y, v.y * v.z,
                           v.z * v.x, v.z * v.y, v.z * v.z);
}

inline double w(const vox::Vector3D& r, const vox::Matrix3x3D& g, double gDet) {
    static const double sigma = 315.0 / (64 * vox::kPiD);
    return sigma * gDet * p((g * r).length());
}

//

AnisotropicPointsToImplicit3::AnisotropicPointsToImplicit3(
                                                           double kernelRadius, double cutOffDensity, double positionSmoothingFactor,
                                                           size_t minNumNeighbors, bool isOutputSdf)
: _kernelRadius(kernelRadius),
_cutOffDensity(cutOffDensity),
_positionSmoothingFactor(positionSmoothingFactor),
_minNumNeighbors(minNumNeighbors),
_isOutputSdf(isOutputSdf) {}

void AnisotropicPointsToImplicit3::convert(
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
    
    LOG(INFO) << "Start converting points to implicit surface.";
    
    const double h = _kernelRadius;
    const double invH = 1 / h;
    const double r = 2.0 * h;
    
    // Mean estimator for cov. mat.
    const auto meanNeighborSearcher =
    vox::PointKdTreeSearcher3::builder().makeShared();
    meanNeighborSearcher->build(points);
    
    LOG(INFO) << "Built neighbor searcher.";
    
    SphSystemData3 meanParticles;
    meanParticles.addParticles(points);
    meanParticles.setNeighborSearcher(meanNeighborSearcher);
    meanParticles.setKernelRadius(r);
    
    // Compute G and xMean
    std::vector<vox::Matrix3x3D> gs(points.size());
    vox::Array1<vox::Vector3D> xMeans(points.size());
    
    vox::parallelFor(vox::kZeroSize, points.size(), [&](size_t i) {
        const auto& x = points[i];
        
        // Compute xMean
        vox::Vector3D xMean;
        double wSum = 0.0;
        size_t numNeighbors = 0;
        const auto getXMean = [&](size_t, const vox::Vector3D& xj) {
            const double wj = wij((x - xj).length(), r);
            wSum += wj;
            xMean += wj * xj;
            ++numNeighbors;
        };
        meanNeighborSearcher->forEachNearbyPoint(x, r, getXMean);
        
        JET_ASSERT(wSum > 0.0);
        xMean /= wSum;
        
        xMeans[i] = lerp(x, xMean, _positionSmoothingFactor);
        
        if (numNeighbors < _minNumNeighbors) {
            const auto g = vox::Matrix3x3D::makeScaleMatrix(invH, invH, invH);
            gs[i] = g;
        } else {
            // Compute covariance matrix
            // We start with small scale matrix (h*h) in order to
            // prevent zero covariance matrix when points are all
            // perfectly lined up.
            auto cov = vox::Matrix3x3D::makeScaleMatrix(h * h, h * h, h * h);
            wSum = 0.0;
            const auto getCov = [&](size_t, const vox::Vector3D& xj) {
                const double wj = wij((xMean - xj).length(), r);
                wSum += wj;
                cov += wj * vvt(xj - xMean);
            };
            meanNeighborSearcher->forEachNearbyPoint(x, r, getCov);
            
            cov /= wSum;
            
            // SVD
            vox::Matrix3x3D u;
            vox::Vector3D v;
            vox::Matrix3x3D w;
            svd(cov, u, v, w);
            
            // Take off the sign
            v.x = std::fabs(v.x);
            v.y = std::fabs(v.y);
            v.z = std::fabs(v.z);
            
            // Constrain Sigma
            const double maxSingularVal = v.max();
            const double kr = 4.0;
            v.x = std::max(v.x, maxSingularVal / kr);
            v.y = std::max(v.y, maxSingularVal / kr);
            v.z = std::max(v.z, maxSingularVal / kr);
            
            const auto invSigma = vox::Matrix3x3D::makeScaleMatrix(1.0 / v);
            
            // Compute G
            const double scale =
            std::pow(v.x * v.y * v.z, 1.0 / 3.0);  // volume preservation
            const vox::Matrix3x3D g = invH * scale * (w * invSigma * u.transposed());
            gs[i] = g;
        }
    });
    
    LOG(INFO) << "Computed G and means.";
    
    // SPH estimator
    meanParticles.setKernelRadius(h);
    meanParticles.updateDensities();
    const auto d = meanParticles.densities();
    const double m = meanParticles.mass();
    
    vox::PointKdTreeSearcher3 meanNeighborSearcher2;
    meanNeighborSearcher2.build(xMeans);
    
    // Compute SDF
    auto temp = output->clone();
    temp->fill([&](const vox::Vector3D& x) {
        double sum = 0.0;
        meanNeighborSearcher2.forEachNearbyPoint(
                                                 x, r,
                                                 [&](size_t i, const vox::Vector3D& neighborPosition) {
            sum += m / d[i] *
            w(neighborPosition - x, gs[i], gs[i].determinant());
        });
        
        return _cutOffDensity - sum;
    }, vox::ExecutionPolicy::kSerial);
    
    LOG(INFO) << "Computed SDF.";
    
    if (_isOutputSdf) {
        FmmLevelSetSolver3 solver;
        solver.reinitialize(*temp, vox::kMaxD, output);
        
        LOG(INFO) << "Completed einitialization.";
    } else {
        temp->swap(output);
    }
    
    LOG(INFO) << "Done converting points to implicit surface.";
}
