//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_HELPER_H_
#define INCLUDE_VDB_HELPER_H_

#include "../src.common/serial.h"
#include "../src.common/parallel.h"
#include "../src.common/array3.h"
#include "../src.common/array_accessor3.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/GridOperators.h>
#include <memory>
#include <vector>

namespace vdb {

template<typename GridType>
class GridHelper{
public:
    using ValueType = typename GridType::ValueType;
    
    typename GridType::Ptr grid;
    
public:
    //! Constructs an empty grid.
    GridHelper(typename GridType::Ptr grid);
    
    template<typename OtherGridType>
    GridHelper<OtherGridType> clone(ValueType init) const;
    
    void swap(GridHelper<GridType> other);
    
    //! Return the transform
    openvdb::math::Transform::Ptr getTransform();
    
    //! Return the const transform
    openvdb::math::Transform::ConstPtr getTransform() const;
    
    openvdb::CoordBBox getBoundingCoordBox();
    
    openvdb::BBoxd getBoudingBox();
    
public:
    //! Returns the grid data at given data point.
    const ValueType& operator()(const openvdb::Coord& pt) const;
    
    //!
    //! \brief Returns the sampled value at given position \p x.
    //!
    //! This function returns the data sampled at arbitrary position \p x.
    //! The sampling function is linear.
    //!
    ValueType sample(const openvdb::Vec3d& x) const;
    
    //!
    //! \brief Returns the sampler function.
    //!
    //! This function returns the data sampler function object. The sampling
    //! function is linear.
    //!
    std::function<ValueType(const openvdb::Vec3d&)> sampler() const;
    
public:
    //! Fills the grid with given value.
    void fill(ValueType value,
              vox::ExecutionPolicy policy = vox::ExecutionPolicy::kParallel);
    
    //! Fills the grid with given position-to-value mapping function.
    void fill(const std::function<ValueType(const openvdb::Vec3d&)>& func,
              vox::ExecutionPolicy policy = vox::ExecutionPolicy::kParallel);
    
    //!
    //! \brief Invokes the given function \p func for each data point.
    //!
    //! This function invokes the given function object \p func for each data
    //! point in serial manner. The input parameters are i and j indices of a
    //! data point. The order of execution is i-first, j-last.
    //!
    void forEachDataPointIndex(
                               const std::function<void(const openvdb::Coord&)>& func) const;
    
    //!
    //! \brief Invokes the given function \p func for each data point
    //! parallelly.
    //!
    //! This function invokes the given function object \p func for each data
    //! point in parallel manner. The input parameters are i and j indices of a
    //! data point. The order of execution can be arbitrary since it's
    //! multi-threaded.
    //!
    void parallelForEachDataPointIndex(
                                       const std::function<void(const openvdb::Coord&)>& func) const;
    
};
//! Shared pointer for the ScalarGrid3 type.

template<typename GridType>
using GridHelperPtr = std::shared_ptr<GridHelper<GridType> >;

template <typename GridType>
void extrapolateToRegion(
                         typename GridType::Ptr input,
                         const vox::ConstArrayAccessor3<char>& valid,
                         unsigned int numberOfIterations,
                         typename GridType::Ptr output);


}  // namespace vox

#include "vdb_helper-inl.h"

#endif  // INCLUDE_VDB_VDB_HELPER_H_
