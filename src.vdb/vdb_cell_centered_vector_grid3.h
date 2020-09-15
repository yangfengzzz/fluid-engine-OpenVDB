//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_CELL_CENTERED_VECTOR_GRID3_H_
#define INCLUDE_VDB_CELL_CENTERED_VECTOR_GRID3_H_

#include "vdb_collocated_vector_grid3.h"
#include <utility>  // just make cpplint happy..
#include <vector>

namespace vdb {

//!
//! \brief 3-D Cell-centered vector grid structure.
//!
//! This class represents 3-D cell-centered vector grid which extends
//! CollocatedVectorGrid3. As its name suggests, the class defines the data
//! point at the center of a grid cell. Thus, the dimension of data points are
//! equal to the dimension of the cells.
//!
class CellCenteredVectorGrid3 final : public CollocatedVectorGrid3 {
public:
    JET_GRID3_TYPE_NAME(CellCenteredVectorGrid3)
    
    class Builder;
    
    //! Constructs zero-sized grid.
    CellCenteredVectorGrid3();
    
    //! Constructs a grid with given resolution, grid spacing, origin and
    //! initial value.
    CellCenteredVectorGrid3(
                            size_t resolutionX,
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
    
    //! Constructs a grid with given resolution, grid spacing, origin and
    //! initial value.
    CellCenteredVectorGrid3(const vox::Size3& resolution,
                            const vox::Vector3D& gridSpacing = vox::Vector3D(1.0, 1.0, 1.0),
                            const vox::Vector3D& origin = vox::Vector3D(),
                            const vox::Vector3D& initialValue = vox::Vector3D());
    
    //! Copy constructor.
    CellCenteredVectorGrid3(const CellCenteredVectorGrid3& other);
    
    //! Returns the actual data point size.
    vox::Size3 dataSize() const override;
    
    //! Returns data position for the grid point at (0, 0, 0).
    //! Note that this is different from origin() since origin() returns
    //! the lower corner point of the bounding box.
    vox::Vector3D dataOrigin() const override;
    
    //! Returns the copy of the grid instance.
    std::shared_ptr<VectorGrid3> clone() const override;
    
    //!
    //! \brief Swaps the contents with the given \p other grid.
    //!
    //! This function swaps the contents of the grid instance with the given
    //! grid object \p other only if \p other has the same type with this grid.
    //!
    void swap(Grid3* other) override;
    
    //! Sets the contents with the given \p other grid.
    void set(const CellCenteredVectorGrid3& other);
    
    //! Sets the contents with the given \p other grid.
    CellCenteredVectorGrid3& operator=(const CellCenteredVectorGrid3& other);
    
    //! Returns builder fox CellCenteredVectorGrid3.
    static Builder builder();
};

//! Shared pointer for the CellCenteredVectorGrid3 type.
typedef std::shared_ptr<CellCenteredVectorGrid3> CellCenteredVectorGrid3Ptr;


//!
//! \brief Front-end to create CellCenteredScalarGrid3 objects step by step.
//!
class CellCenteredVectorGrid3::Builder final : public VectorGridBuilder3 {
public:
    //! Returns builder with resolution.
    Builder& withResolution(const vox::Size3& resolution);
    
    //! Returns builder with resolution.
    Builder& withResolution(
                            size_t resolutionX,
                            size_t resolutionY,
                            size_t resolutionZ);
    
    //! Returns builder with grid spacing.
    Builder& withGridSpacing(const vox::Vector3D& gridSpacing);
    
    //! Returns builder with grid spacing.
    Builder& withGridSpacing(
                             double gridSpacingX, double gridSpacingY, double gridSpacingZ);
    
    //! Returns builder with grid origin.
    Builder& withOrigin(const vox::Vector3D& gridOrigin);
    
    //! Returns builder with grid origin.
    Builder& withOrigin(
                        double gridOriginX, double gridOriginY, double gridOriginZ);
    
    //! Returns builder with initial value.
    Builder& withInitialValue(const vox::Vector3D& initialVal);
    
    //! Returns builder with initial value.
    Builder& withInitialValue(
                              double initialValX, double initialValY, double initialValZ);
    
    //! Builds CellCenteredScalarGrid3 instance.
    CellCenteredVectorGrid3 build() const;
    
    //! Builds shared pointer of CellCenteredVectorGrid3 instance.
    CellCenteredVectorGrid3Ptr makeShared() const;
    
    //!
    //! \brief Builds shared pointer of CellCenteredVectorGrid3 instance.
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

#endif  // INCLUDE_VDB_CELL_CENTERED_VECTOR_GRID3_H_
