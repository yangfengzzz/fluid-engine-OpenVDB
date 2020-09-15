//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_COLLOCATED_VECTOR_GRID3_H_
#define INCLUDE_VDB_COLLOCATED_VECTOR_GRID3_H_

#include "vdb_samplers3.h"
#include "vdb_vector_grid3.h"
#include <vector>

namespace vdb {

//! \brief Abstract base class for 3-D collocated vector grid structure.
class CollocatedVectorGrid3 : public VectorGrid3 {
public:
    //! Constructs an empty grid.
    CollocatedVectorGrid3();
    
    //! Default destructor.
    virtual ~CollocatedVectorGrid3();
    
    const openvdb::Vec3dGrid::Ptr getGrid() const {return _grid;}
    
    //! Returns the actual data point size.
    virtual vox::Size3 dataSize() const = 0;
    
    //!
    //! \brief Returns data position for the grid point at (0, 0, 0).
    //!
    //! Note that this is different from origin() since origin() returns
    //! the lower corner point of the bounding box.
    //!
    virtual vox::Vector3D dataOrigin() const = 0;
    
    //! Returns the grid data at given data point.
    openvdb::Vec3d operator()(const openvdb::Coord& coord);
    
    //! Returns divergence at data point location.
    double divergenceAtDataPoint(const openvdb::Coord& coord) const;
    
    //! Returns curl at data point location.
    vox::Vector3D curlAtDataPoint(const openvdb::Coord& coord) const;
    
    //! Returns the read-write data array accessor.
    openvdb::Vec3dGrid::Accessor dataAccessor();
    
    //! Returns the read-only data array accessor.
    openvdb::Vec3dGrid::ConstAccessor constDataAccessor() const;
    
    //! Returns the function that maps data point to its position.
    DataPositionFunc dataPosition() const;
    
    //!
    //! \brief Invokes the given function \p func for each data point.
    //!
    //! This function invokes the given function object \p func for each data
    //! point in serial manner. The input parameters are i and j indices of a
    //! data point. The order of execution is i-first, j-last.
    //!
    void forEachDataPointIndex(
                               const std::function<void(const openvdb::Coord& coord)>& func) const;
    
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
                                       const std::function<void(const openvdb::Coord& coord)>& func) const;
    
    //! Fills the grid with given value.
    void fill(const vox::Vector3D& value,
              vox::ExecutionPolicy policy = vox::ExecutionPolicy::kParallel) override;
    
    //! Fills the grid with given function.
    void fill(const std::function<vox::Vector3D(const vox::Vector3D&)>& func,
              vox::ExecutionPolicy policy = vox::ExecutionPolicy::kParallel) override;
    
    // VectorField3 implementations
    
    //! Returns sampled value at given position \p x.
    vox::Vector3D sample(const vox::Vector3D& x) const override;
    
    //! Returns divergence at given position \p x.
    double divergence(const vox::Vector3D& x) const override;
    
    //! Returns curl at given position \p x.
    vox::Vector3D curl(const vox::Vector3D& x) const override;
    
    //!
    //! \brief Returns the sampler function.
    //!
    //! This function returns the data sampler function object. The sampling
    //! function is linear.
    //!
    std::function<vox::Vector3D(const vox::Vector3D&)> sampler() const override;
    
protected:
    //! Swaps the data storage and predefined samplers with given grid.
    void swapCollocatedVectorGrid(CollocatedVectorGrid3* other);
    
    //! Sets the data storage and predefined samplers with given grid.
    void setCollocatedVectorGrid(const CollocatedVectorGrid3& other);
    
private:
    openvdb::Vec3dGrid::Ptr _grid;
    LinearGridSampler<openvdb::Vec3dGrid> _linearSampler;
    std::function<vox::Vector3D(const vox::Vector3D&)> _sampler;
    
    void onResize(
                  const vox::Size3& resolution,
                  const vox::Vector3D& gridSpacing,
                  const vox::Vector3D& origin,
                  const vox::Vector3D& initialValue) final;
    
    void resetSampler();
};

//! Shared pointer for the CollocatedVectorGrid3 type.
typedef std::shared_ptr<CollocatedVectorGrid3> CollocatedVectorGrid3Ptr;

}  // namespace vox

#endif  // INCLUDE_VDB_COLLOCATED_VECTOR_GRID3_H_
