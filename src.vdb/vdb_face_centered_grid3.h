//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_FACE_CENTERED_GRID3_H_
#define INCLUDE_VDB_FACE_CENTERED_GRID3_H_

#include "vdb_samplers3.h"
#include "vdb_vector_grid3.h"
#include <memory>
#include <utility>  // just make cpplint happy..
#include <vector>

namespace vdb {

//!
//! \brief 3-D face-centered (a.k.a MAC or staggered) grid.
//!
//! This class implements face-centered grid which is also known as
//! marker-and-cell (MAC) or staggered grid. This vector grid stores each vector
//! component at face center. Thus, u, v, and w components are not collocated.
//!
class FaceCenteredGrid3 final : public VectorGrid3 {
public:
    JET_GRID3_TYPE_NAME(FaceCenteredGrid3)
    
    class Builder;
    
    //! Constructs empty grid.
    FaceCenteredGrid3();
    
    //! Resizes the grid using given parameters.
    FaceCenteredGrid3(size_t resolutionX,
                      size_t resolutionY,
                      size_t resolutionZ,
                      double gridSpacingX = 1.0,
                      double gridSpacingY = 1.0,
                      double gridSpacingZ = 1.0,
                      double originX = 0.0,
                      double originY = 0.0,
                      double originZ = 0.0,
                      double initialValueU = 0.0,
                      double initialValueV = 0.0,
                      double initialValueW = 0.0);
    
    //! Resizes the grid using given parameters.
    FaceCenteredGrid3(const vox::Size3& resolution,
                      const vox::Vector3D& gridSpacing = vox::Vector3D(1.0, 1.0, 1.0),
                      const vox::Vector3D& origin = vox::Vector3D(),
                      const vox::Vector3D& initialValue = vox::Vector3D());
    
    //! Copy constructor.
    FaceCenteredGrid3(const FaceCenteredGrid3& other);
    
    const openvdb::DoubleGrid::Ptr getUGrid() const {return _dataU;}
    
    const openvdb::DoubleGrid::Ptr getVGrid() const {return _dataV;}
    
    const openvdb::DoubleGrid::Ptr getWGrid() const {return _dataW;}
    
    //!
    //! \brief Swaps the contents with the given \p other grid.
    //!
    //! This function swaps the contents of the grid instance with the given
    //! grid object \p other only if \p other has the same type with this grid.
    //!
    void swap(Grid3* other) override;
    
    //! Sets the contents with the given \p other grid.
    void set(const FaceCenteredGrid3& other);
    
    //! Sets the contents with the given \p other grid.
    FaceCenteredGrid3& operator=(const FaceCenteredGrid3& other);
    
    //! Returns u-value at given data point.
    double u(const openvdb::Coord& coord) const;
    
    //! Returns v-value at given data point.
    double v(const openvdb::Coord& coord) const;
    
    //! Returns w-value at given data point.
    double w(const openvdb::Coord& coord) const;
    
    //! Returns interpolated value at cell center.
    vox::Vector3D valueAtCellCenter(int i, int j, int k) const;
    
    //! Returns divergence at cell-center location.
    double divergenceAtCellCenter(int i, int j, int k) const;
    
    //! Returns curl at cell-center location.
    vox::Vector3D curlAtCellCenter(int i, int j, int k) const;
    
    //! Returns u data accessor.
    openvdb::DoubleGrid::Accessor uAccessor();
    
    //! Returns read-only u data accessor.
    openvdb::DoubleGrid::ConstAccessor uConstAccessor() const;
    
    //! Returns v data accessor.
    openvdb::DoubleGrid::Accessor vAccessor();
    
    //! Returns read-only v data accessor.
    openvdb::DoubleGrid::ConstAccessor vConstAccessor() const;
    
    //! Returns w data accessor.
    openvdb::DoubleGrid::Accessor wAccessor();
    
    //! Returns read-only w data accessor.
    openvdb::DoubleGrid::ConstAccessor wConstAccessor() const;
    
    //! Returns function object that maps u data point to its actual position.
    DataPositionFunc uPosition() const;
    
    //! Returns function object that maps v data point to its actual position.
    DataPositionFunc vPosition() const;
    
    //! Returns function object that maps w data point to its actual position.
    DataPositionFunc wPosition() const;
    
    //! Returns data size of the u component.
    vox::Size3 uSize() const;
    
    //! Returns data size of the v component.
    vox::Size3 vSize() const;
    
    //! Returns data size of the w component.
    vox::Size3 wSize() const;
    
    //!
    //! \brief Returns u-data position for the grid point at (0, 0, 0).
    //!
    //! Note that this is different from origin() since origin() returns
    //! the lower corner point of the bounding box.
    //!
    vox::Vector3D uOrigin() const;
    
    //!
    //! \brief Returns v-data position for the grid point at (0, 0, 0).
    //!
    //! Note that this is different from origin() since origin() returns
    //! the lower corner point of the bounding box.
    //!
    vox::Vector3D vOrigin() const;
    
    //!
    //! \brief Returns w-data position for the grid point at (0, 0, 0).
    //!
    //! Note that this is different from origin() since origin() returns
    //! the lower corner point of the bounding box.
    //!
    vox::Vector3D wOrigin() const;
    
    //! Fills the grid with given value.
    void fill(const vox::Vector3D& value,
              vox::ExecutionPolicy policy = vox::ExecutionPolicy::kParallel) override;
    
    //! Fills the grid with given function.
    void fill(const std::function<vox::Vector3D(const vox::Vector3D&)>& func,
              vox::ExecutionPolicy policy = vox::ExecutionPolicy::kParallel) override;
    
    //! Returns the copy of the grid instance.
    std::shared_ptr<VectorGrid3> clone() const override;
    
    //!
    //! \brief Invokes the given function \p func for each u-data point.
    //!
    //! This function invokes the given function object \p func for each u-data
    //! point in serial manner. The input parameters are i and j indices of a
    //! u-data point. The order of execution is i-first, j-last.
    //!
    void forEachUIndex(
                       const std::function<void(const openvdb::Coord& coord)>& func) const;
    
    //!
    //! \brief Invokes the given function \p func for each u-data point
    //! parallelly.
    //!
    //! This function invokes the given function object \p func for each u-data
    //! point in parallel manner. The input parameters are i and j indices of a
    //! u-data point. The order of execution can be arbitrary since it's
    //! multi-threaded.
    //!
    void parallelForEachUIndex(
                               const std::function<void(const openvdb::Coord& coord)>& func) const;
    
    //!
    //! \brief Invokes the given function \p func for each v-data point.
    //!
    //! This function invokes the given function object \p func for each v-data
    //! point in serial manner. The input parameters are i and j indices of a
    //! v-data point. The order of execution is i-first, j-last.
    //!
    void forEachVIndex(
                       const std::function<void(const openvdb::Coord& coord)>& func) const;
    
    //!
    //! \brief Invokes the given function \p func for each v-data point
    //! parallelly.
    //!
    //! This function invokes the given function object \p func for each v-data
    //! point in parallel manner. The input parameters are i and j indices of a
    //! v-data point. The order of execution can be arbitrary since it's
    //! multi-threaded.
    //!
    void parallelForEachVIndex(
                               const std::function<void(const openvdb::Coord& coord)>& func) const;
    
    //!
    //! \brief Invokes the given function \p func for each w-data point.
    //!
    //! This function invokes the given function object \p func for each w-data
    //! point in serial manner. The input parameters are i and j indices of a
    //! w-data point. The order of execution is i-first, j-last.
    //!
    void forEachWIndex(
                       const std::function<void(const openvdb::Coord& coord)>& func) const;
    
    //!
    //! \brief Invokes the given function \p func for each w-data point
    //! parallelly.
    //!
    //! This function invokes the given function object \p func for each w-data
    //! point in parallel manner. The input parameters are i and j indices of a
    //! w-data point. The order of execution can be arbitrary since it's
    //! multi-threaded.
    //!
    void parallelForEachWIndex(
                               const std::function<void(const openvdb::Coord& coord)>& func) const;
    
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
    
    //! Returns builder fox FaceCenteredGrid3.
    static Builder builder();
    
protected:
    // VectorGrid3 implementations
    void onResize(const vox::Size3& resolution,
                  const vox::Vector3D& gridSpacing,
                  const vox::Vector3D& origin,
                  const vox::Vector3D& initialValue) final;
    
private:
    openvdb::DoubleGrid::Ptr _dataU;
    openvdb::DoubleGrid::Ptr _dataV;
    openvdb::DoubleGrid::Ptr _dataW;
    vox::Vector3D _dataOriginU = vox::Vector3D();
    vox::Vector3D _dataOriginV = vox::Vector3D();
    vox::Vector3D _dataOriginW = vox::Vector3D();
    
    LinearGridSampler<openvdb::DoubleGrid> _uLinearSampler;
    LinearGridSampler<openvdb::DoubleGrid> _vLinearSampler;
    LinearGridSampler<openvdb::DoubleGrid> _wLinearSampler;
    std::function<vox::Vector3D(const vox::Vector3D&)> _sampler;
    
    void resetSampler();
};

//! Shared pointer type for the FaceCenteredGrid3.
typedef std::shared_ptr<FaceCenteredGrid3> FaceCenteredGrid3Ptr;

//!
//! \brief Front-end to create CellCenteredScalarGrid3 objects step by step.
//!
class FaceCenteredGrid3::Builder final : public VectorGridBuilder3 {
public:
    //! Returns builder with resolution.
    Builder& withResolution(const vox::Size3& resolution);
    
    //! Returns builder with resolution.
    Builder& withResolution(size_t resolutionX,
                            size_t resolutionY,
                            size_t resolutionZ);
    
    //! Returns builder with grid spacing.
    Builder& withGridSpacing(const vox::Vector3D& gridSpacing);
    
    //! Returns builder with grid spacing.
    Builder& withGridSpacing(double gridSpacingX,
                             double gridSpacingY,
                             double gridSpacingZ);
    
    //! Returns builder with grid origin.
    Builder& withOrigin(const vox::Vector3D& gridOrigin);
    
    //! Returns builder with grid origin.
    Builder& withOrigin(double gridOriginX,
                        double gridOriginY,
                        double gridOriginZ);
    
    //! Returns builder with initial value.
    Builder& withInitialValue(const vox::Vector3D& initialVal);
    
    //! Returns builder with initial value.
    Builder& withInitialValue(double initialValX,
                              double initialValY,
                              double initialValZ);
    
    //! Builds CellCenteredScalarGrid3 instance.
    FaceCenteredGrid3 build() const;
    
    //! Builds shared pointer of FaceCenteredGrid3 instance.
    FaceCenteredGrid3Ptr makeShared() const;
    
    //!
    //! \brief Builds shared pointer of FaceCenteredGrid3 instance.
    //!
    //! This is an overriding function that implements VectorGridBuilder3.
    //!
    VectorGrid3Ptr build(const vox::Size3& resolution,
                         const vox::Vector3D& gridSpacing,
                         const vox::Vector3D& gridOrigin,
                         const vox::Vector3D& initialVal) const override;
    
private:
    vox::Size3 _resolution{1, 1, 1};
    vox::Vector3D _gridSpacing{1, 1, 1};
    vox::Vector3D _gridOrigin{0, 0, 0};
    vox::Vector3D _initialVal{0, 0, 0};
};

}  // namespace vox

#endif  // INCLUDE_VDB_FACE_CENTERED_GRID3_H_
