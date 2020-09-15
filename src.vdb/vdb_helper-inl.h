//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "../src.common/pch.h"

#include "../src.common/parallel.h"
#include "../src.common/fdm_utils.h"
#include <openvdb/tools/Interpolation.h>

using namespace vdb;

template<typename GridType>
GridHelper<GridType>::GridHelper(typename GridType::Ptr grid):
grid(grid){
}

template<typename GridType>
template<typename OtherGridType>
GridHelper<OtherGridType> GridHelper<GridType>::clone(ValueType init) const{
    GridHelper<OtherGridType> copy_grid(OtherGridType::create(init));
    copy_grid.grid->setTransform(getTransform()->copy() );
    copy_grid.grid->topologyUnion(*grid);
    return copy_grid;
}

template<typename GridType>
void GridHelper<GridType>::swap(GridHelper<GridType> other){
    std::swap(grid, other.grid);
}

template<typename GridType>
openvdb::math::Transform::Ptr GridHelper<GridType>::getTransform(){
    return grid->transformPtr();
}

template<typename GridType>
openvdb::math::Transform::ConstPtr GridHelper<GridType>::getTransform() const{
    return grid->transformPtr();
}

template<typename GridType>
openvdb::CoordBBox GridHelper<GridType>::getBoundingCoordBox(){
    return grid->evalActiveVoxelBoundingBox();
}

template<typename GridType>
openvdb::BBoxd GridHelper<GridType>::getBoudingBox(){
    return grid->transform().indexToWorld(getBoundingCoordBox());
}
//------------------------------------------------------------------------------
template<typename GridType>
const typename GridHelper<GridType>::ValueType&
GridHelper<GridType>::operator()(const openvdb::Coord& pt) const {
    return grid->tree().getValue(pt);
}

template<typename GridType>
typename GridHelper<GridType>::ValueType
GridHelper<GridType>::sample(const openvdb::Vec3d& x) const {
    openvdb::tools::GridSampler<GridType, openvdb::tools::BoxSampler> sampler(*grid);
    return sampler.wsSample(x);
}

template<typename GridType>
std::function<typename GridHelper<GridType>::ValueType(const openvdb::Vec3d&)>
GridHelper<GridType>::sampler() const {
    return [&](const openvdb::Vec3d& x)->GridHelper<GridType>::ValueType{
        openvdb::tools::GridSampler<GridType, openvdb::tools::BoxSampler> sampler(*grid);
        return sampler.wsSample(x);
    };
}
//------------------------------------------------------------------------------
template<typename GridType>
void GridHelper<GridType>::fill(GridHelper<GridType>::ValueType value,
                                vox::ExecutionPolicy policy) {
    using IterRange = openvdb::tree::IteratorRange<typename GridType::ValueOnIter>;
    IterRange range(grid->beginValueOn());
    tbb::parallel_for(range, [&](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            typename GridType::ValueOnIter iter = r.iterator();
            iter.setValue(value);
        }
    });
}

template<typename GridType>
void GridHelper<GridType>::fill(const std::function<GridHelper<GridType>::ValueType(const openvdb::Vec3d&)>& func,
                                vox::ExecutionPolicy policy) {
    using IterRange = openvdb::tree::IteratorRange<typename GridType::ValueOnIter>;
    IterRange range(grid->beginValueOn());
    tbb::parallel_for(range, [&](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            typename GridType::ValueOnIter iter = r.iterator();
            const openvdb::Vec3d& xyz = grid->transform().indexToWorld(iter.getCoord());
            iter.setValue(func(xyz));
        }
    });
}

template<typename GridType>
void GridHelper<GridType>::forEachDataPointIndex(
                                                 const std::function<void(const openvdb::Coord&)>& func) const {
    for (typename GridType::ValueOnIter iter = grid->beginValueOn(); iter; ++iter) {
        func(iter.getCoord());
    }
}

template<typename GridType>
void GridHelper<GridType>::parallelForEachDataPointIndex(
                                                         const std::function<void(const openvdb::Coord&)>& func) const {
    using IterRange = openvdb::tree::IteratorRange<typename GridType::ValueOnIter>;
    IterRange range(grid->beginValueOn());
    tbb::parallel_for(range, [&func](IterRange& r) {
        // Iterate over a subrange of the leaf iterator's iteration space.
        for ( ; r; ++r) {
            typename GridType::ValueOnIter iter = r.iterator();
            func(iter.getCoord());
        }
    });
}
//----------------------------------------------------------------
template <typename GridType>
void vdb::extrapolateToRegion(typename GridType::Ptr input,
                              const vox::ConstArrayAccessor3<char>& valid,
                              unsigned int numberOfIterations,
                              typename GridType::Ptr output){
    const vox::Size3 size = valid.size();
    
    vox::Array3<char> valid0(size);
    vox::Array3<char> valid1(size);
    
    valid0.forEachIndex([&](uint i, uint j, uint k) {
        valid0(i, j, k) = valid(i, j, k);
        openvdb::Coord coord(i, j, k);
        output->tree().setValueOnly(coord, input->tree().getValue(coord));
    });
    
    for (unsigned int iter = 0; iter < numberOfIterations; ++iter) {
        valid0.forEachIndex([&](uint i, uint j, uint k) {
            typename GridType::ValueType sum = openvdb::zeroVal<typename GridType::ValueType>();
            unsigned int count = 0;
            
            if (!valid0(i, j, k)) {
                if (i + 1 < size.x && valid0(i + 1, j, k)) {
                    sum += output->tree().getValue(openvdb::Coord(i + 1, j, k));
                    ++count;
                }
                
                if (i > 0 && valid0(i - 1, j, k)) {
                    sum += output->tree().getValue(openvdb::Coord(i - 1, j, k));
                    ++count;
                }
                
                if (j + 1 < size.y && valid0(i, j + 1, k)) {
                    sum += output->tree().getValue(openvdb::Coord(i, j + 1, k));
                    ++count;
                }
                
                if (j > 0 && valid0(i, j - 1, k)) {
                    sum += output->tree().getValue(openvdb::Coord(i, j - 1, k));
                    ++count;
                }
                
                if (k + 1 < size.z && valid0(i, j, k + 1)) {
                    sum += output->tree().getValue(openvdb::Coord(i, j, k + 1));
                    ++count;
                }
                
                if (k > 0 && valid0(i, j, k - 1)) {
                    sum += output->tree().getValue(openvdb::Coord(i, j, k - 1));
                    ++count;
                }
                
                if (count > 0) {
                    output->tree().setValueOn(openvdb::Coord(i, j, k),
                                              sum / static_cast<double>(count));
                    valid1(i, j, k) = 1;
                }
            } else {
                valid1(i, j, k) = 1;
            }
        });
        
        valid0.swap(valid1);
    }
}
