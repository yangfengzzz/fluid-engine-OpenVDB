//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_ARRAY_SAMPLERS3_H_
#define INCLUDE_VDB_ARRAY_SAMPLERS3_H_

#include "../src.common/vector3.h"
#include "../src.common/size3.h"
#include <openvdb/openvdb.h>
#include <functional>

namespace vdb {

//!
//! \brief 3-D nearest grid sampler class.
//!
//! This class provides nearest sampling interface for a given 3-D array.
//!
//! \tparam GridType - The value type to sample.
//!
template <typename GridType>
class NearestGridSampler final {
public:
    //!
    //! \brief      Constructs a sampler using array accessor, spacing between
    //!     the elements, and the position of the first array element.
    //!
    //! \param[in]  gridSpacing The grid spacing.
    //! \param[in]  gridOrigin  The grid origin.
    //!
    explicit NearestGridSampler(
                                const typename GridType::Ptr& grid,
                                const vox::Size3& resolution,
                                const vox::Vector3D& gridSpacing,
                                const vox::Vector3D& gridOrigin);
    
    //! Copy constructor.
    NearestGridSampler(const NearestGridSampler& other);
    
    //! Returns sampled value at point \p pt.
    typename GridType::ValueType operator()(const vox::Vector3D& pt) const;
    
    //! Returns the nearest array index for point \p x.
    void getCoordinate(const vox::Vector3D& pt, openvdb::Coord* index) const;
    
    //! Returns a funtion object that wraps this instance.
    std::function<typename GridType::ValueType(const vox::Vector3D&)> functor() const;
    
private:
    vox::Size3 _resolution;
    vox::Vector3D _gridSpacing;
    vox::Vector3D _origin;
    typename GridType::Ptr _grid;
};

//!
//! \brief 3-D linear grid sampler class.
//!
//! This class provides linear sampling interface for a given 2-D array.
//!
//! \tparam GridType - The value type to sample.
//!
template <typename GridType>
class LinearGridSampler final {
public:
    //!
    //! \brief      Constructs a sampler using array accessor, spacing between
    //!     the elements, and the position of the first array element.
    //!
    //! \param[in]  grid    The array accessor.
    //!
    explicit LinearGridSampler(
                               const typename GridType::Ptr& grid,
                               const vox::Size3& resolution,
                               const vox::Vector3D& gridSpacing,
                               const vox::Vector3D& gridOrigin);
    
    //! Copy constructor.
    LinearGridSampler(const LinearGridSampler& other);
    
    //! Returns sampled value at point \p pt.
    typename GridType::ValueType operator()(const vox::Vector3D& pt) const;
    
    //! Returns the indices of points and their sampling weight for given point.
    void getCoordinatesAndWeights(
                                  const vox::Vector3D& pt,
                                  std::array<openvdb::Coord, 8>* indices,
                                  std::array<double, 8>* weights) const;
    
    //! Returns the indices of points and their gradient of sampling weight for
    //! given point.
    void getCoordinatesAndGradientWeights(
                                          const vox::Vector3D& pt,
                                          std::array<openvdb::Coord, 8>* indices,
                                          std::array<vox::Vector3D, 8>* weights) const;
    
    //! Returns a funtion object that wraps this instance.
    std::function<typename GridType::ValueType(const vox::Vector3D&)> functor() const;
    
private:
    vox::Size3 _resolution;
    vox::Vector3D _gridSpacing;
    vox::Vector3D _invGridSpacing;
    vox::Vector3D _origin;
    typename GridType::Ptr _grid;
};

//!
//! \brief 3-D cubic grid sampler class.
//!
//! This class provides cubic sampling interface for a given 3-D array.
//!
//! \tparam GridType - The grid type to sample.
//!
template <typename GridType>
class CubicGridSampler final {
public:
    //!
    //! \brief      Constructs a sampler using array accessor, spacing between
    //!     the elements, and the position of the first array element.
    //!
    //! \param[in]  gridSpacing The grid spacing.
    //! \param[in]  gridOrigin  The grid origin.
    //!
    explicit CubicGridSampler(
                              const typename GridType::Ptr& grid,
                              const vox::Size3& resolution,
                              const vox::Vector3D& gridSpacing,
                              const vox::Vector3D& gridOrigin);
    
    //! Copy constructor.
    CubicGridSampler(const CubicGridSampler& other);
    
    //! Returns sampled value at point \p pt.
    typename GridType::ValueType operator()(const vox::Vector3D& pt) const;
    
    //! Returns a funtion object that wraps this instance.
    std::function<typename GridType::ValueType(const vox::Vector3D&)> functor() const;
    
private:
    vox::Size3 _resolution;
    vox::Vector3D _gridSpacing;
    vox::Vector3D _origin;
    typename GridType::Ptr _grid;
};

}  // namespace vox

#include "vdb_samplers3-inl.h"

#endif  // INCLUDE_VDB_ARRAY_SAMPLERS3_H_
