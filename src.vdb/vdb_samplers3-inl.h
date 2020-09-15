//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_DETAIL_ARRAY_SAMPLERS3_INL_H_
#define INCLUDE_VDB_DETAIL_ARRAY_SAMPLERS3_INL_H_

#include "../src.common/pch.h"
#include "../src.common/macros.h"
#include "vdb_math_utils.h"

#include <openvdb/tools/Interpolation.h>
#include <algorithm>
#include <functional>
#include <limits>

namespace vdb {

template <typename GridType>
NearestGridSampler<GridType>::NearestGridSampler(
                                                 const typename GridType::Ptr& grid,
                                                 const vox::Size3& resolution,
                                                 const vox::Vector3D& gridSpacing,
                                                 const vox::Vector3D& gridOrigin):
_grid(grid){
    _resolution = resolution;
    _gridSpacing = gridSpacing;
    _origin = gridOrigin;
}

template <typename GridType>
NearestGridSampler<GridType>::NearestGridSampler(
                                                 const NearestGridSampler& other):
_grid(other._grid){
    _resolution = other._resolution;
    _gridSpacing = other._gridSpacing;
    _origin = other._origin;
}

template <typename GridType>
typename GridType::ValueType NearestGridSampler<GridType>::operator()(const vox::Vector3D& x) const {
    ssize_t i, j, k;
    double fx, fy, fz;
    
    JET_ASSERT(_gridSpacing.x > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.y > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.z > std::numeric_limits<double>::epsilon());
    vox::Vector3D normalizedX = (x - _origin) / _gridSpacing;
    
    ssize_t iSize = static_cast<ssize_t>(_resolution.x);
    ssize_t jSize = static_cast<ssize_t>(_resolution.y);
    ssize_t kSize = static_cast<ssize_t>(_resolution.z);
    
    vox::getBarycentric(normalizedX.x, 0, iSize - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, jSize - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, kSize - 1, &k, &fz);
    
    i = std::min(static_cast<ssize_t>(i + fx + 0.5), iSize - 1);
    j = std::min(static_cast<ssize_t>(j + fy + 0.5), jSize - 1);
    k = std::min(static_cast<ssize_t>(k + fz + 0.5), kSize - 1);
    
    return _grid->tree().getValue(openvdb::Coord(i, j, k));
}

template <typename GridType>
void NearestGridSampler<GridType>::getCoordinate(
                                                 const vox::Vector3D& x,
                                                 openvdb::Coord* index) const {
    ssize_t i, j, k;
    double fx, fy, fz;
    
    JET_ASSERT(_gridSpacing.x > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.y > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.z > std::numeric_limits<double>::epsilon());
    vox::Vector3D normalizedX = (x - _origin) / _gridSpacing;
    
    ssize_t iSize = static_cast<ssize_t>(_resolution.x);
    ssize_t jSize = static_cast<ssize_t>(_resolution.y);
    ssize_t kSize = static_cast<ssize_t>(_resolution.z);
    
    vox::getBarycentric(normalizedX.x, 0, iSize - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, jSize - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, kSize - 1, &k, &fz);
    
    index->x() = std::min(static_cast<ssize_t>(i + fx + 0.5), iSize - 1);
    index->y() = std::min(static_cast<ssize_t>(j + fy + 0.5), jSize - 1);
    index->z() = std::min(static_cast<ssize_t>(k + fz + 0.5), kSize - 1);
}

template <typename GridType>
std::function<typename GridType::ValueType(const vox::Vector3D&)>

NearestGridSampler<GridType>::functor() const {
    NearestGridSampler sampler(*this);
    return std::bind(
                     &NearestGridSampler::operator(), sampler, std::placeholders::_1);
}
//----------------------------------------------------------------------------

template <typename GridType>
LinearGridSampler<GridType>::LinearGridSampler(
                                               const typename GridType::Ptr& grid,
                                               const vox::Size3& resolution,
                                               const vox::Vector3D& gridSpacing,
                                               const vox::Vector3D& gridOrigin)
{
    _resolution = resolution;
    _gridSpacing = gridSpacing;
    _invGridSpacing = 1.0 / _gridSpacing;
    _origin = gridOrigin;
    _grid = grid;
}

template <typename GridType>
LinearGridSampler<GridType>::LinearGridSampler(
                                               const LinearGridSampler& other){
    _resolution = other._resolution;
    _gridSpacing = other._gridSpacing;
    _invGridSpacing = other._invGridSpacing;
    _origin = other._origin;
    _grid = other._grid;
}

template <typename GridType>
typename GridType::ValueType LinearGridSampler<GridType>::operator()(const vox::Vector3D& x) const {
    ssize_t i, j, k;
    double fx, fy, fz;
    
    JET_ASSERT(_gridSpacing.x > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.y > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.z > std::numeric_limits<double>::epsilon());
    vox::Vector3D normalizedX = (x - _origin) / _gridSpacing;
    
    ssize_t iSize = static_cast<ssize_t>(_resolution.x);
    ssize_t jSize = static_cast<ssize_t>(_resolution.y);
    ssize_t kSize = static_cast<ssize_t>(_resolution.z);
    
    vox::getBarycentric(normalizedX.x, 0, iSize - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, jSize - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, kSize - 1, &k, &fz);
    
    ssize_t ip1 = std::min(i + 1, iSize - 1);
    ssize_t jp1 = std::min(j + 1, jSize - 1);
    ssize_t kp1 = std::min(k + 1, kSize - 1);
    
    return vox::trilerp(
                        _grid->tree().getValue(openvdb::Coord(i, j, k)),
                        _grid->tree().getValue(openvdb::Coord(ip1, j, k)),
                        _grid->tree().getValue(openvdb::Coord(i, jp1, k)),
                        _grid->tree().getValue(openvdb::Coord(ip1, jp1, k)),
                        _grid->tree().getValue(openvdb::Coord(i, j, kp1)),
                        _grid->tree().getValue(openvdb::Coord(ip1, j, kp1)),
                        _grid->tree().getValue(openvdb::Coord(i, jp1, kp1)),
                        _grid->tree().getValue(openvdb::Coord(ip1, jp1, kp1)),
                        fx,
                        fy,
                        fz);
}

template <typename GridType>
void LinearGridSampler<GridType>::getCoordinatesAndWeights(
                                                           const vox::Vector3D& x,
                                                           std::array<openvdb::Coord, 8>* indices,
                                                           std::array<double, 8>* weights) const {
    ssize_t i, j, k;
    double fx, fy, fz;
    
    JET_ASSERT(
               _gridSpacing.x > 0.0 && _gridSpacing.y > 0.0 && _gridSpacing.z > 0.0);
    
    const vox::Vector3D normalizedX = (x - _origin) * _invGridSpacing;
    
    ssize_t iSize = static_cast<ssize_t>(_resolution.x);
    ssize_t jSize = static_cast<ssize_t>(_resolution.y);
    ssize_t kSize = static_cast<ssize_t>(_resolution.z);
    
    vox::getBarycentric(normalizedX.x, 0, iSize - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, jSize - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, kSize - 1, &k, &fz);
    
    ssize_t ip1 = std::min(i + 1, iSize - 1);
    ssize_t jp1 = std::min(j + 1, jSize - 1);
    ssize_t kp1 = std::min(k + 1, kSize - 1);
    
    (*indices)[0] = openvdb::Coord(i, j, k);
    (*indices)[1] = openvdb::Coord(ip1, j, k);
    (*indices)[2] = openvdb::Coord(i, jp1, k);
    (*indices)[3] = openvdb::Coord(ip1, jp1, k);
    (*indices)[4] = openvdb::Coord(i, j, kp1);
    (*indices)[5] = openvdb::Coord(ip1, j, kp1);
    (*indices)[6] = openvdb::Coord(i, jp1, kp1);
    (*indices)[7] = openvdb::Coord(ip1, jp1, kp1);
    
    (*weights)[0] = (1 - fx) * (1 - fy) * (1 - fz);
    (*weights)[1] = fx * (1 - fy) * (1 - fz);
    (*weights)[2] = (1 - fx) * fy * (1 - fz);
    (*weights)[3] = fx * fy * (1 - fz);
    (*weights)[4] = (1 - fx) * (1 - fy) * fz;
    (*weights)[5] = fx * (1 - fy) * fz;
    (*weights)[6] = (1 - fx) * fy * fz;
    (*weights)[7] = fx * fy * fz;
}

template <typename GridType>
void LinearGridSampler<GridType>::getCoordinatesAndGradientWeights(
                                                                   const vox::Vector3D& x,
                                                                   std::array<openvdb::Coord, 8>* indices,
                                                                   std::array<vox::Vector3D, 8>* weights) const{
    ssize_t i, j, k;
    double fx, fy, fz;
    
    JET_ASSERT(
               _gridSpacing.x > 0.0 && _gridSpacing.y > 0.0 && _gridSpacing.z > 0.0);
    
    const vox::Vector3D normalizedX = (x - _origin) * _invGridSpacing;
    
    ssize_t iSize = static_cast<ssize_t>(_resolution.x);
    ssize_t jSize = static_cast<ssize_t>(_resolution.y);
    ssize_t kSize = static_cast<ssize_t>(_resolution.z);
    
    vox::getBarycentric(normalizedX.x, 0, iSize - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, jSize - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, kSize - 1, &k, &fz);
    
    ssize_t ip1 = std::min(i + 1, iSize - 1);
    ssize_t jp1 = std::min(j + 1, jSize - 1);
    ssize_t kp1 = std::min(k + 1, kSize - 1);
    
    (*indices)[0] = openvdb::Coord(i, j, k);
    (*indices)[1] = openvdb::Coord(ip1, j, k);
    (*indices)[2] = openvdb::Coord(i, jp1, k);
    (*indices)[3] = openvdb::Coord(ip1, jp1, k);
    (*indices)[4] = openvdb::Coord(i, j, kp1);
    (*indices)[5] = openvdb::Coord(ip1, j, kp1);
    (*indices)[6] = openvdb::Coord(i, jp1, kp1);
    (*indices)[7] = openvdb::Coord(ip1, jp1, kp1);
    
    (*weights)[0] = vox::Vector3D(
                                  -_invGridSpacing.x * (1 - fy) * (1 - fz),
                                  -_invGridSpacing.y * (1 - fx) * (1 - fz),
                                  -_invGridSpacing.z * (1 - fx) * (1 - fy));
    (*weights)[1] = vox::Vector3D(
                                  _invGridSpacing.x * (1 - fy) * (1 - fz),
                                  fx * (-_invGridSpacing.y) * (1 - fz),
                                  fx * (1 - fy) * (-_invGridSpacing.z));
    (*weights)[2] = vox::Vector3D(
                                  (-_invGridSpacing.x) * fy * (1 - fz),
                                  (1 - fx) * _invGridSpacing.y * (1 - fz),
                                  (1 - fx) * fy * (-_invGridSpacing.z));
    (*weights)[3] = vox::Vector3D(
                                  _invGridSpacing.x * fy * (1 - fz),
                                  fx * _invGridSpacing.y * (1 - fz),
                                  fx * fy * (-_invGridSpacing.z));
    (*weights)[4] = vox::Vector3D(
                                  (-_invGridSpacing.x) * (1 - fy) * fz,
                                  (1 - fx) * (-_invGridSpacing.y) * fz,
                                  (1 - fx) * (1 - fy) * _invGridSpacing.z);
    (*weights)[5] = vox::Vector3D(
                                  _invGridSpacing.x * (1 - fy) * fz,
                                  fx * (-_invGridSpacing.y) * fz,
                                  fx * (1 - fy) * _invGridSpacing.z);
    (*weights)[6] = vox::Vector3D(
                                  (-_invGridSpacing.x) * fy * fz,
                                  (1 - fx) * _invGridSpacing.y * fz,
                                  (1 - fx) * fy * _invGridSpacing.z);
    (*weights)[7] = vox::Vector3D(
                                  _invGridSpacing.x * fy * fz,
                                  fx * _invGridSpacing.y * fz,
                                  fx * fy * _invGridSpacing.z);
}

template <typename GridType>
std::function<typename GridType::ValueType(const vox::Vector3D&)>
LinearGridSampler<GridType>::functor() const {
    LinearGridSampler sampler(*this);
    return std::bind(
                     &LinearGridSampler::operator(), sampler, std::placeholders::_1);
}

//-----------------------------------------------------------------------------
template <typename GridType>
CubicGridSampler<GridType>::CubicGridSampler(
                                             const typename GridType::Ptr& grid,
                                             const vox::Size3& resolution,
                                             const vox::Vector3D& gridSpacing,
                                             const vox::Vector3D& gridOrigin) :
_grid(grid){
    _resolution = resolution;
    _gridSpacing = gridSpacing;
    _origin = gridOrigin;
}


template <typename GridType>
CubicGridSampler<GridType>::CubicGridSampler(
                                             const CubicGridSampler& other):
_grid(other._grid){
    _resolution = other._resolution;
    _gridSpacing = other._gridSpacing;
    _origin = other._origin;
}

template <typename GridType>
typename GridType::ValueType CubicGridSampler<GridType>::operator()(const vox::Vector3D& x) const {
    ssize_t i, j, k;
    ssize_t iSize = static_cast<ssize_t>(_resolution.x);
    ssize_t jSize = static_cast<ssize_t>(_resolution.y);
    ssize_t kSize = static_cast<ssize_t>(_resolution.z);
    double fx, fy, fz;
    
    JET_ASSERT(_gridSpacing.x > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.y > std::numeric_limits<double>::epsilon() &&
               _gridSpacing.z > std::numeric_limits<double>::epsilon());
    vox::Vector3D normalizedX = (x - _origin) / _gridSpacing;
    
    vox::getBarycentric(normalizedX.x, 0, iSize - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, jSize - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, kSize - 1, &k, &fz);
    
    ssize_t is[4] = {
        std::max(i - 1, vox::kZeroSSize),
        i,
        std::min(i + 1, iSize - 1),
        std::min(i + 2, iSize - 1)
    };
    ssize_t js[4] = {
        std::max(j - 1, vox::kZeroSSize),
        j,
        std::min(j + 1, jSize - 1),
        std::min(j + 2, jSize - 1)
    };
    ssize_t ks[4] = {
        std::max(k - 1, vox::kZeroSSize),
        k,
        std::min(k + 1, kSize - 1),
        std::min(k + 2, kSize - 1)
    };
    
    typename GridType::ValueType kValues[4];
    
    for (int kk = 0; kk < 4; ++kk) {
        typename GridType::ValueType jValues[4];
        
        for (int jj = 0; jj < 4; ++jj) {
            jValues[jj] = vox::monotonicCatmullRom(
                                                   _grid->tree().getValue(openvdb::Coord(is[0], js[jj], ks[kk])),
                                                   _grid->tree().getValue(openvdb::Coord(is[1], js[jj], ks[kk])),
                                                   _grid->tree().getValue(openvdb::Coord(is[2], js[jj], ks[kk])),
                                                   _grid->tree().getValue(openvdb::Coord(is[3], js[jj], ks[kk])),
                                                   fx);
        }
        
        kValues[kk] = vox::monotonicCatmullRom(
                                               jValues[0], jValues[1], jValues[2], jValues[3], fy);
    }
    
    return vox::monotonicCatmullRom(
                                    kValues[0], kValues[1], kValues[2], kValues[3], fz);
}

template <typename GridType>
std::function<typename GridType::ValueType(const vox::Vector3D&)> CubicGridSampler<GridType>::functor() const {
    CubicGridSampler sampler(*this);
    return std::bind(
                     &CubicGridSampler::operator(), sampler, std::placeholders::_1);
}

}  // namespace vox

#endif  // INCLUDE_VDB_DETAIL_ARRAY_SAMPLERS3_INL_H_
