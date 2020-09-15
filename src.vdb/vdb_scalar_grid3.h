//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_SCALAR_GRID3_H_
#define INCLUDE_VDB_SCALAR_GRID3_H_

#include "vdb_grid3.h"
#include "../src.common/scalar_field3.h"
#include "../src.common/array_accessor3.h"
#include "vdb_samplers3.h"
#include <memory>
#include <vector>

namespace vdb {

//! Abstract base class for 3-D scalar grid structure.
class ScalarGrid3 : public vox::ScalarField3, public Grid3 {
public:
    //! Constructs an empty grid.
    ScalarGrid3();
    
    //! Default destructor.
    virtual ~ScalarGrid3();
    
    const openvdb::DoubleGrid::Ptr getGrid() const {return _grid;}
    
    //!
    //! \brief Returns the size of the grid data.
    //!
    //! This function returns the size of the grid data which is not necessarily
    //! equal to the grid resolution if the data is not stored at cell-center.
    //!
    virtual vox::Size3 dataSize() const = 0;
    
    //!
    //! \brief Returns the origin of the grid data.
    //!
    //! This function returns data position for the grid point at (0, 0, 0).
    //! Note that this is different from origin() since origin() returns
    //! the lower corner point of the bounding box.
    //!
    virtual vox::Vector3D dataOrigin() const = 0;
    
    //! Returns the copy of the grid instance.
    virtual std::shared_ptr<ScalarGrid3> clone() const = 0;
    
    //! Clears the contents of the grid.
    void clear();
    
    //! Resizes the grid using given parameters.
    void resize(
                size_t resolutionX,
                size_t resolutionY,
                size_t resolutionZ,
                double gridSpacingX = 1.0,
                double gridSpacingY = 1.0,
                double gridSpacingZ = 1.0,
                double originX = 0.0,
                double originY = 0.0,
                double originZ = 0.0,
                double initialValue = 0.0);
    
    //! Resizes the grid using given parameters.
    void resize(
                const vox::Size3& resolution,
                const vox::Vector3D& gridSpacing = vox::Vector3D(1, 1, 1),
                const vox::Vector3D& origin = vox::Vector3D(),
                double initialValue = 0.0);
    
    //! Resizes the grid using given parameters.
    void resize(
                double gridSpacingX,
                double gridSpacingY,
                double gridSpacingZ,
                double originX,
                double originY,
                double originZ);
    
    //! Resizes the grid using given parameters.
    void resize(const vox::Vector3D& gridSpacing,
                const vox::Vector3D& origin);
    
    //! Returns the grid data at given data point.
    double operator()(const openvdb::Coord& coord) const;
    
    //! Returns the gradient vector at given data point.
    vox::Vector3D gradientAtDataPoint(const openvdb::Coord& coord) const;
    
    //! Returns the Laplacian at given data point.
    double laplacianAtDataPoint(const openvdb::Coord& coord) const;
    
    //! Returns the read-write data array accessor.
    openvdb::DoubleGrid::Accessor dataAccessor();
    
    //! Returns the read-only data array accessor.
    openvdb::DoubleGrid::ConstAccessor constDataAccessor() const;
    
    //! Returns the function that maps data point to its position.
    DataPositionFunc dataPosition() const;
    
    //! Fills the grid with given value.
    void fill(double value,
              vox::ExecutionPolicy policy = vox::ExecutionPolicy::kParallel);
    
    //! Fills the grid with given position-to-value mapping function.
    void fill(const std::function<double(const vox::Vector3D&)>& func,
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
    
    // ScalarField3 implementations
    
    //!
    //! \brief Returns the sampled value at given position \p x.
    //!
    //! This function returns the data sampled at arbitrary position \p x.
    //! The sampling function is linear.
    //!
    double sample(const vox::Vector3D& x) const override;
    
    //!
    //! \brief Returns the sampler function.
    //!
    //! This function returns the data sampler function object. The sampling
    //! function is linear.
    //!
    std::function<double(const vox::Vector3D&)> sampler() const override;
    
    //! Returns the gradient vector at given position \p x.
    vox::Vector3D gradient(const vox::Vector3D& x) const override;
    
    //! Returns the Laplacian at given position \p x.
    double laplacian(const vox::Vector3D& x) const override;
    
    //! Fetches the data into a continuous linear array.
    void getData(vox::ArrayAccessor3<double> data) const;
    
    //! Sets the data from a continuous linear array.
    void setData(const vox::ConstArrayAccessor3<double> data);
    
protected:
    //! Swaps the data storage and predefined samplers with given grid.
    void swapScalarGrid(ScalarGrid3* other);
    
    //! Sets the data storage and predefined samplers with given grid.
    void setScalarGrid(const ScalarGrid3& other);
    
private:
    openvdb::DoubleGrid::Ptr _grid;
    LinearGridSampler<openvdb::DoubleGrid> _linearSampler;
    std::function<double(const vox::Vector3D&)> _sampler;
    
    void resetSampler();
};

//! Shared pointer for the ScalarGrid3 type.
typedef std::shared_ptr<ScalarGrid3> ScalarGrid3Ptr;

//! Abstract base class for 3-D scalar grid builder.
class ScalarGridBuilder3 {
public:
    //! Creates a builder.
    ScalarGridBuilder3();
    
    //! Default destructor.
    virtual ~ScalarGridBuilder3();
    
    //! Returns 3-D scalar grid with given parameters.
    virtual ScalarGrid3Ptr build(
                                 const vox::Size3& resolution,
                                 const vox::Vector3D& gridSpacing,
                                 const vox::Vector3D& gridOrigin,
                                 double initialVal) const = 0;
};

//! Shared pointer for the ScalarGridBuilder3 type.
typedef std::shared_ptr<ScalarGridBuilder3> ScalarGridBuilder3Ptr;

}  // namespace vox

#endif  // INCLUDE_VDB_SCALAR_GRID3_H_
