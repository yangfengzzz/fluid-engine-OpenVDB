//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "../src.common/pch.h"
#include "vdb_samplers3.h"
#include "vdb_face_centered_grid3.h"
#include "../src.common/parallel.h"
#include "../src.common/serial.h"
#include <openvdb/tools/Interpolation.h>
#include <algorithm>
#include <utility>  // just make cpplint happy..
#include <vector>

using namespace vdb;

FaceCenteredGrid3::FaceCenteredGrid3()
: _dataOriginU(0.0, 0.5, 0.5),
_dataOriginV(0.5, 0.0, 0.5),
_dataOriginW(0.5, 0.5, 0.0),
_uLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataU, uSize(),
                                                       vox::Vector3D(1,1,1),
                                                       _dataOriginU)),
_vLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataV, vSize(),
                                                       vox::Vector3D(1,1,1),
                                                       _dataOriginV)),
_wLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataW, wSize(),
                                                       vox::Vector3D(1,1,1),
                                                       _dataOriginW)){}

FaceCenteredGrid3::FaceCenteredGrid3(size_t resolutionX,
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
                                     double initialValueW)
: FaceCenteredGrid3(vox::Size3(resolutionX,
                               resolutionY,
                               resolutionZ),
                    vox::Vector3D(gridSpacingX,
                                  gridSpacingY,
                                  gridSpacingZ),
                    vox::Vector3D(originX,
                                  originY,
                                  originZ),
                    vox::Vector3D(initialValueU,
                                  initialValueV,
                                  initialValueW)) {
}

FaceCenteredGrid3::FaceCenteredGrid3(const vox::Size3& resolution,
                                     const vox::Vector3D& gridSpacing,
                                     const vox::Vector3D& origin,
                                     const vox::Vector3D& initialValue)
:_uLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataU, uSize(),
                                                        vox::Vector3D(1,1,1),
                                                        _dataOriginU)),
_vLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataV, vSize(),
                                                       vox::Vector3D(1,1,1),
                                                       _dataOriginV)),
_wLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataW, wSize(),
                                                       vox::Vector3D(1,1,1),
                                                       _dataOriginW))
{
    resize(resolution, gridSpacing, origin, initialValue);
}

FaceCenteredGrid3::FaceCenteredGrid3(const FaceCenteredGrid3& other)
:_uLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataU, uSize(),
                                                        vox::Vector3D(1,1,1),
                                                        _dataOriginU)),
_vLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataV, vSize(),
                                                       vox::Vector3D(1,1,1),
                                                       _dataOriginV)),
_wLinearSampler(LinearGridSampler<openvdb::DoubleGrid>(_dataW, wSize(),
                                                       vox::Vector3D(1,1,1),
                                                       _dataOriginW)) {
    set(other);
}

void FaceCenteredGrid3::swap(Grid3* other) {
    FaceCenteredGrid3* sameType = dynamic_cast<FaceCenteredGrid3*>(other);
    
    if (sameType != nullptr) {
        swapGrid(sameType);
        
        _dataU.swap(sameType->_dataU);
        _dataV.swap(sameType->_dataV);
        _dataW.swap(sameType->_dataW);
        std::swap(_dataOriginU, sameType->_dataOriginU);
        std::swap(_dataOriginV, sameType->_dataOriginV);
        std::swap(_dataOriginW, sameType->_dataOriginW);
        std::swap(_uLinearSampler, sameType->_uLinearSampler);
        std::swap(_vLinearSampler, sameType->_vLinearSampler);
        std::swap(_wLinearSampler, sameType->_wLinearSampler);
        std::swap(_sampler, sameType->_sampler);
    }
}

void FaceCenteredGrid3::set(const FaceCenteredGrid3& other) {
    setGrid(other);
    
    _dataU = other._dataU->deepCopy();
    _dataV = other._dataV->deepCopy();
    _dataW = other._dataW->deepCopy();
    _dataOriginU = other._dataOriginU;
    _dataOriginV = other._dataOriginV;
    _dataOriginW = other._dataOriginW;
    
    resetSampler();
}

FaceCenteredGrid3& FaceCenteredGrid3::operator=(
                                                const FaceCenteredGrid3& other)
{
    set(other);
    return *this;
}

double FaceCenteredGrid3::u(const openvdb::Coord& coord) const{
    return _dataU->tree().getValue(coord);
}

double FaceCenteredGrid3::v(const openvdb::Coord& coord) const{
    return _dataV->tree().getValue(coord);
}

double FaceCenteredGrid3::w(const openvdb::Coord& coord) const{
    return _dataW->tree().getValue(coord);
}

vox::Vector3D FaceCenteredGrid3::valueAtCellCenter(int i,
                                                   int j,
                                                   int k) const {
    JET_ASSERT(i < resolution().x && j < resolution().y && k < resolution().z);
    
    return 0.5 * vox::Vector3D(_dataU->tree().getValue(openvdb::Coord(i, j, k))
                               + _dataU->tree().getValue(openvdb::Coord(i + 1, j, k)),
                               _dataV->tree().getValue(openvdb::Coord(i, j, k))
                               + _dataV->tree().getValue(openvdb::Coord(i, j + 1, k)),
                               _dataW->tree().getValue(openvdb::Coord(i, j, k))
                               + _dataW->tree().getValue(openvdb::Coord(i, j, k + 1)));
}

double FaceCenteredGrid3::divergenceAtCellCenter(int i,
                                                 int j,
                                                 int k) const {
    JET_ASSERT(i < resolution().x && j < resolution().y && k < resolution().z);
    
    const vox::Vector3D& gs = gridSpacing();
    
    double leftU = _dataU->tree().getValue(openvdb::Coord(i, j, k));
    double rightU = _dataU->tree().getValue(openvdb::Coord(i + 1, j, k));
    double bottomV = _dataV->tree().getValue(openvdb::Coord(i, j, k));
    double topV = _dataV->tree().getValue(openvdb::Coord(i, j + 1, k));
    double backW = _dataW->tree().getValue(openvdb::Coord(i, j, k));
    double frontW = _dataW->tree().getValue(openvdb::Coord(i, j, k + 1));
    
    return (rightU - leftU) / gs.x
    + (topV - bottomV) / gs.y
    + (frontW - backW) / gs.z;
}

vox::Vector3D FaceCenteredGrid3::curlAtCellCenter(int i,
                                                  int j,
                                                  int k) const {
    const vox::Size3& res = resolution();
    const vox::Vector3D& gs = gridSpacing();
    
    JET_ASSERT(i < res.x && j < res.y && k < res.z);
    
    vox::Vector3D left = valueAtCellCenter((i > 0) ? i - 1 : i, j, k);
    vox::Vector3D right = valueAtCellCenter((i + 1 < res.x) ? i + 1 : i, j, k);
    vox::Vector3D down = valueAtCellCenter(i, (j > 0) ? j - 1 : j, k);
    vox::Vector3D up = valueAtCellCenter(i, (j + 1 < res.y) ? j + 1 : j, k);
    vox::Vector3D back = valueAtCellCenter(i, j, (k > 0) ? k - 1 : k);
    vox::Vector3D front = valueAtCellCenter(i, j, (k + 1 < res.z) ? k + 1 : k);
    
    double Fx_ym = down.x;
    double Fx_yp = up.x;
    double Fx_zm = back.x;
    double Fx_zp = front.x;
    
    double Fy_xm = left.y;
    double Fy_xp = right.y;
    double Fy_zm = back.y;
    double Fy_zp = front.y;
    
    double Fz_xm = left.z;
    double Fz_xp = right.z;
    double Fz_ym = down.z;
    double Fz_yp = up.z;
    
    return vox::Vector3D(
                         0.5 * (Fz_yp - Fz_ym) / gs.y
                         - 0.5 * (Fy_zp - Fy_zm) / gs.z,
                         0.5 * (Fx_zp - Fx_zm) / gs.z
                         - 0.5 * (Fz_xp - Fz_xm) / gs.x,
                         0.5 * (Fy_xp - Fy_xm) / gs.x
                         - 0.5 * (Fx_yp - Fx_ym) / gs.y);
}


openvdb::DoubleGrid::Accessor FaceCenteredGrid3::uAccessor(){
    return _dataU->getAccessor();
}

openvdb::DoubleGrid::ConstAccessor FaceCenteredGrid3::uConstAccessor() const{
    return _dataU->getConstAccessor();
}

openvdb::DoubleGrid::Accessor FaceCenteredGrid3::vAccessor(){
    return _dataV->getAccessor();
}

openvdb::DoubleGrid::ConstAccessor FaceCenteredGrid3::vConstAccessor() const{
    return _dataV->getConstAccessor();
}

openvdb::DoubleGrid::Accessor FaceCenteredGrid3::wAccessor(){
    return _dataW->getAccessor();
}

openvdb::DoubleGrid::ConstAccessor FaceCenteredGrid3::wConstAccessor() const{
    return _dataW->getConstAccessor();
}

VectorGrid3::DataPositionFunc FaceCenteredGrid3::uPosition() const {
    vox::Vector3D h = gridSpacing();
    
    return [this, h](const openvdb::Coord& coord) -> vox::Vector3D {
        return _dataOriginU + h * vox::Vector3D(coord.x(),
                                                coord.y(),
                                                coord.z());
    };
}

VectorGrid3::DataPositionFunc FaceCenteredGrid3::vPosition() const {
    vox::Vector3D h = gridSpacing();
    
    return [this, h](const openvdb::Coord& coord) -> vox::Vector3D {
        return _dataOriginV + h * vox::Vector3D(coord.x(),
                                                coord.y(),
                                                coord.z());
    };
}

VectorGrid3::DataPositionFunc FaceCenteredGrid3::wPosition() const {
    vox::Vector3D h = gridSpacing();
    
    return [this, h](const openvdb::Coord& coord) -> vox::Vector3D {
        return _dataOriginW + h * vox::Vector3D(coord.x(),
                                                coord.y(),
                                                coord.z());
    };
}

vox::Size3 FaceCenteredGrid3::uSize() const {
    if (resolution() != vox::Size3(0, 0, 0)) {
        return resolution() + vox::Size3(1, 0, 0);
    } else {
        return vox::Size3(0, 0, 0);
    }
}

vox::Size3 FaceCenteredGrid3::vSize() const {
    if (resolution() != vox::Size3(0, 0, 0)) {
        return resolution() + vox::Size3(0, 1, 0);
    } else {
        return vox::Size3(0, 0, 0);
    }
}

vox::Size3 FaceCenteredGrid3::wSize() const {
    if (resolution() != vox::Size3(0, 0, 0)) {
        return resolution() + vox::Size3(0, 0, 1);
    } else {
        return vox::Size3(0, 0, 0);
    }
}

vox::Vector3D FaceCenteredGrid3::uOrigin() const { return _dataOriginU; }

vox::Vector3D FaceCenteredGrid3::vOrigin() const { return _dataOriginV; }

vox::Vector3D FaceCenteredGrid3::wOrigin() const { return _dataOriginW; }

void FaceCenteredGrid3::fill(const vox::Vector3D& value,
                             vox::ExecutionPolicy policy) {
    vox::parallelFor(
                     vox::kZeroSize, uSize().x,
                     vox::kZeroSize, uSize().y,
                     vox::kZeroSize, uSize().z,
                     [this, value](uint i, uint j, uint k) {
        _dataU->tree().setValueOnly(openvdb::Coord(i, j, k), value.x);
    }, policy);
    
    vox::parallelFor(
                     vox::kZeroSize, vSize().x,
                     vox::kZeroSize, vSize().y,
                     vox::kZeroSize, vSize().z,
                     [this, value](uint i, uint j, uint k) {
        _dataV->tree().setValueOnly(openvdb::Coord(i, j, k), value.y);
    }, policy);
    
    vox::parallelFor(
                     vox::kZeroSize, wSize().x,
                     vox::kZeroSize, wSize().y,
                     vox::kZeroSize, wSize().z,
                     [this, value](uint i, uint j, uint k) {
        _dataW->tree().setValueOnly(openvdb::Coord(i, j, k), value.z);
    }, policy);
}

void FaceCenteredGrid3::fill(
                             const std::function<vox::Vector3D(const vox::Vector3D&)>& func,
                             vox::ExecutionPolicy policy) {
    DataPositionFunc uPos = uPosition();
    vox::parallelFor(vox::kZeroSize, uSize().x,
                     vox::kZeroSize, uSize().y,
                     vox::kZeroSize, uSize().z,
                     [this, &func, &uPos](uint i, uint j, uint k) {
        _dataU->tree().setValueOnly(openvdb::Coord(i, j, k),
                                    func(uPos(openvdb::Coord(i, j, k))).x);
    }, policy);
    
    DataPositionFunc vPos = vPosition();
    vox::parallelFor(vox::kZeroSize, vSize().x,
                     vox::kZeroSize, vSize().y,
                     vox::kZeroSize, vSize().z,
                     [this, &func, &vPos](uint i, uint j, uint k) {
        _dataV->tree().setValueOnly(openvdb::Coord(i, j, k),
                                    func(vPos(openvdb::Coord(i, j, k))).y);
    }, policy);
    
    DataPositionFunc wPos = wPosition();
    vox::parallelFor(vox::kZeroSize, wSize().x,
                     vox::kZeroSize, wSize().y,
                     vox::kZeroSize, wSize().z,
                     [this, &func, &wPos](uint i, uint j, uint k) {
        _dataW->tree().setValueOnly(openvdb::Coord(i, j, k),
                                    func(wPos(openvdb::Coord(i, j, k))).z);
    }, policy);
}

std::shared_ptr<VectorGrid3> FaceCenteredGrid3::clone() const {
    return CLONE_W_CUSTOM_DELETER(FaceCenteredGrid3);
}

void FaceCenteredGrid3::forEachUIndex(
                                      const std::function<void(const openvdb::Coord& coord)>& func) const {
    for (openvdb::DoubleGrid::ValueOnIter iter = _dataU->beginValueOn(); iter; ++iter) {
        func(iter.getCoord());
    }
}

void FaceCenteredGrid3::parallelForEachUIndex(
                                              const std::function<void(const openvdb::Coord& coord)>& func) const {
    using IterRange = openvdb::tree::IteratorRange<openvdb::DoubleGrid::ValueOnIter>;
    IterRange range(_dataU->beginValueOn());
    tbb::parallel_for(range, [&func](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            openvdb::DoubleGrid::ValueOnIter iter = r.iterator();
            func(iter.getCoord());
        }
    });
}

void FaceCenteredGrid3::forEachVIndex(
                                      const std::function<void(const openvdb::Coord& coord)>& func) const {
    for (openvdb::DoubleGrid::ValueOnIter iter = _dataV->beginValueOn(); iter; ++iter) {
        func(iter.getCoord());
    }
}

void FaceCenteredGrid3::parallelForEachVIndex(
                                              const std::function<void(const openvdb::Coord& coord)>& func) const {
    using IterRange = openvdb::tree::IteratorRange<openvdb::DoubleGrid::ValueOnIter>;
    IterRange range(_dataV->beginValueOn());
    tbb::parallel_for(range, [&func](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            openvdb::DoubleGrid::ValueOnIter iter = r.iterator();
            func(iter.getCoord());
        }
    });
}

void FaceCenteredGrid3::forEachWIndex(
                                      const std::function<void(const openvdb::Coord& coord)>& func) const {
    for (openvdb::DoubleGrid::ValueOnIter iter = _dataW->beginValueOn(); iter; ++iter) {
        func(iter.getCoord());
    }
}

void FaceCenteredGrid3::parallelForEachWIndex(
                                              const std::function<void(const openvdb::Coord& coord)>& func) const {
    using IterRange = openvdb::tree::IteratorRange<openvdb::DoubleGrid::ValueOnIter>;
    IterRange range(_dataW->beginValueOn());
    tbb::parallel_for(range, [&func](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            openvdb::DoubleGrid::ValueOnIter iter = r.iterator();
            func(iter.getCoord());
        }
    });
}

vox::Vector3D FaceCenteredGrid3::sample(const vox::Vector3D& x) const {
    return _sampler(x);
}

std::function<vox::Vector3D(const vox::Vector3D&)> FaceCenteredGrid3::sampler() const {
    return _sampler;
}

double FaceCenteredGrid3::divergence(const vox::Vector3D& x) const {
    vox::Size3 res = resolution();
    ssize_t i, j, k;
    double fx, fy, fz;
    vox::Vector3D cellCenterOrigin = origin() + 0.5 * gridSpacing();
    
    vox::Vector3D normalizedX = (x - cellCenterOrigin) / gridSpacing();
    
    vox::getBarycentric(normalizedX.x, 0, static_cast<ssize_t>(res.x) - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, static_cast<ssize_t>(res.y) - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, static_cast<ssize_t>(res.z) - 1, &k, &fz);
    
    std::array<openvdb::Coord, 8> indices;
    std::array<double, 8> weights;
    
    indices[0] = openvdb::Coord(i, j, k);
    indices[1] = openvdb::Coord(i + 1, j, k);
    indices[2] = openvdb::Coord(i, j + 1, k);
    indices[3] = openvdb::Coord(i + 1, j + 1, k);
    indices[4] = openvdb::Coord(i, j, k + 1);
    indices[5] = openvdb::Coord(i + 1, j, k + 1);
    indices[6] = openvdb::Coord(i, j + 1, k + 1);
    indices[7] = openvdb::Coord(i + 1, j + 1, k + 1);
    
    weights[0] = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
    weights[1] = fx * (1.0 - fy) * (1.0 - fz);
    weights[2] = (1.0 - fx) * fy * (1.0 - fz);
    weights[3] = fx * fy * (1.0 - fz);
    weights[4] = (1.0 - fx) * (1.0 - fy) * fz;
    weights[5] = fx * (1.0 - fy) * fz;
    weights[6] = (1.0 - fx) * fy * fz;
    weights[7] = fx * fy * fz;
    
    double result = 0.0;
    
    for (int n = 0; n < 8; ++n) {
        result += weights[n] * divergenceAtCellCenter(
                                                      indices[n].x(), indices[n].y(), indices[n].z());
    }
    
    return result;
}

vox::Vector3D FaceCenteredGrid3::curl(const vox::Vector3D& x) const {
    vox::Size3 res = resolution();
    ssize_t i, j, k;
    double fx, fy, fz;
    vox::Vector3D cellCenterOrigin = origin() + 0.5 * gridSpacing();
    
    vox::Vector3D normalizedX = (x - cellCenterOrigin) / gridSpacing();
    
    vox::getBarycentric(normalizedX.x, 0, static_cast<ssize_t>(res.x) - 1, &i, &fx);
    vox::getBarycentric(normalizedX.y, 0, static_cast<ssize_t>(res.y) - 1, &j, &fy);
    vox::getBarycentric(normalizedX.z, 0, static_cast<ssize_t>(res.z) - 1, &k, &fz);
    
    std::array<openvdb::Coord, 8> indices;
    std::array<double, 8> weights;
    
    indices[0] = openvdb::Coord(i, j, k);
    indices[1] = openvdb::Coord(i + 1, j, k);
    indices[2] = openvdb::Coord(i, j + 1, k);
    indices[3] = openvdb::Coord(i + 1, j + 1, k);
    indices[4] = openvdb::Coord(i, j, k + 1);
    indices[5] = openvdb::Coord(i + 1, j, k + 1);
    indices[6] = openvdb::Coord(i, j + 1, k + 1);
    indices[7] = openvdb::Coord(i + 1, j + 1, k + 1);
    
    weights[0] = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
    weights[1] = fx * (1.0 - fy) * (1.0 - fz);
    weights[2] = (1.0 - fx) * fy * (1.0 - fz);
    weights[3] = fx * fy * (1.0 - fz);
    weights[4] = (1.0 - fx) * (1.0 - fy) * fz;
    weights[5] = fx * (1.0 - fy) * fz;
    weights[6] = (1.0 - fx) * fy * fz;
    weights[7] = fx * fy * fz;
    
    vox::Vector3D result;
    
    for (int n = 0; n < 8; ++n) {
        result += weights[n] *
        curlAtCellCenter(indices[n].x(), indices[n].y(), indices[n].z());
    }
    
    return result;
}

void FaceCenteredGrid3::onResize(const vox::Size3& resolution,
                                 const vox::Vector3D& gridSpacing,
                                 const vox::Vector3D& origin,
                                 const vox::Vector3D& initialValue) {
    _dataU = openvdb::DoubleGrid::create(initialValue.x);
    _dataV = openvdb::DoubleGrid::create(initialValue.y);
    _dataW = openvdb::DoubleGrid::create(initialValue.z);
    
    _dataOriginU = origin + 0.5 * vox::Vector3D(0.0, gridSpacing.y,
                                                gridSpacing.z);
    _dataOriginV = origin + 0.5 * vox::Vector3D(gridSpacing.x, 0.0,
                                                gridSpacing.z);
    _dataOriginW = origin + 0.5 * vox::Vector3D(gridSpacing.x,
                                                gridSpacing.y, 0.0);
    
    resetSampler();
}

void FaceCenteredGrid3::resetSampler() {
    LinearGridSampler<openvdb::DoubleGrid> uSampler(_dataU, uSize(),
                                                    gridSpacing(), _dataOriginU);
    LinearGridSampler<openvdb::DoubleGrid> vSampler(_dataV, vSize(),
                                                    gridSpacing(), _dataOriginV);
    LinearGridSampler<openvdb::DoubleGrid> wSampler(_dataW, wSize(),
                                                    gridSpacing(), _dataOriginW);
    
    _uLinearSampler = uSampler;
    _vLinearSampler = vSampler;
    _wLinearSampler = wSampler;
    
    _sampler = [uSampler, vSampler, wSampler](const vox::Vector3D& x) -> vox::Vector3D {
        double u = uSampler(x);
        double v = vSampler(x);
        double w = wSampler(x);
        return vox::Vector3D(u, v, w);
    };
}

//-------------------------------------------------------------------------
FaceCenteredGrid3::Builder FaceCenteredGrid3::builder() { return Builder(); }

FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withResolution(
                                                                       const vox::Size3& resolution) {
    _resolution = resolution;
    return *this;
}

FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withResolution(
                                                                       size_t resolutionX, size_t resolutionY, size_t resolutionZ) {
    _resolution.x = resolutionX;
    _resolution.y = resolutionY;
    _resolution.z = resolutionZ;
    return *this;
}

FaceCenteredGrid3::Builder&
FaceCenteredGrid3::Builder::withGridSpacing(
                                            const vox::Vector3D& gridSpacing) {
    _gridSpacing = gridSpacing;
    return *this;
}

FaceCenteredGrid3::Builder&
FaceCenteredGrid3::Builder::withGridSpacing(
                                            double gridSpacingX,
                                            double gridSpacingY,
                                            double gridSpacingZ) {
    _gridSpacing.x = gridSpacingX;
    _gridSpacing.y = gridSpacingY;
    _gridSpacing.z = gridSpacingZ;
    return *this;
}

FaceCenteredGrid3::Builder&
FaceCenteredGrid3::Builder::withOrigin(
                                       const vox::Vector3D& gridOrigin) {
    _gridOrigin = gridOrigin;
    return *this;
}

FaceCenteredGrid3::Builder&
FaceCenteredGrid3::Builder::withOrigin(
                                       double gridOriginX,
                                       double gridOriginY,
                                       double gridOriginZ) {
    _gridOrigin.x = gridOriginX;
    _gridOrigin.y = gridOriginY;
    _gridOrigin.z = gridOriginZ;
    return *this;
}

FaceCenteredGrid3::Builder&
FaceCenteredGrid3::Builder::withInitialValue(
                                             const vox::Vector3D& initialVal) {
    _initialVal = initialVal;
    return *this;
}

FaceCenteredGrid3::Builder&
FaceCenteredGrid3::Builder::withInitialValue(
                                             double initialValX,
                                             double initialValY,
                                             double initialValZ) {
    _initialVal.x = initialValX;
    _initialVal.y = initialValY;
    _initialVal.z = initialValZ;
    return *this;
}

FaceCenteredGrid3 FaceCenteredGrid3::Builder::build() const {
    return FaceCenteredGrid3(_resolution,
                             _gridSpacing,
                             _gridOrigin,
                             _initialVal);
}

FaceCenteredGrid3Ptr
FaceCenteredGrid3::Builder::makeShared() const {
    return std::shared_ptr<FaceCenteredGrid3>(
                                              new FaceCenteredGrid3(_resolution,
                                                                    _gridSpacing,
                                                                    _gridOrigin,
                                                                    _initialVal),
                                              [](FaceCenteredGrid3* obj) { delete obj; });
}

VectorGrid3Ptr
FaceCenteredGrid3::Builder::build(const vox::Size3& resolution,
                                  const vox::Vector3D& gridSpacing,
                                  const vox::Vector3D& gridOrigin,
                                  const vox::Vector3D& initialVal) const {
    return std::shared_ptr<FaceCenteredGrid3>(
                                              new FaceCenteredGrid3(resolution,
                                                                    gridSpacing,
                                                                    gridOrigin,
                                                                    initialVal),
                                              [](FaceCenteredGrid3* obj) { delete obj; });
}
