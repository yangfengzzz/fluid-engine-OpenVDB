//  Copyright © 2020 Feng Yang. All rights reserved.
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_VDB_GRID3_H_
#define INCLUDE_VDB_GRID3_H_

#include "../src.common/bounding_box3.h"
#include "../src.common/size3.h"
#include "../src.common/parallel.h"

#include <openvdb/openvdb.h>
#include <functional>
#include <memory>
#include <string>
#include <utility>  // just make cpplint happy..
#include <vector>

namespace vdb {

//!
//! \brief Abstract base class for 3-D cartesian grid structure.
//!
//! This class represents 3-D cartesian grid structure. This class is an
//! abstract base class and does not store any data. The class only stores the
//! shape of the grid. The grid structure is axis-aligned and can have different
//! grid spacing per axis.
//!
class Grid3{
public:
    //! Function type for mapping data index to actual position.
    typedef std::function<vox::Vector3D(openvdb::Coord)> DataPositionFunc;
    
    //! Constructs an empty grid.
    Grid3();
    
    //! Default destructor.
    virtual ~Grid3();
    
    //! Returns the type name of derived grid.
    virtual std::string typeName() const = 0;
    
    //! Returns the grid resolution.
    const vox::Size3& resolution() const;
    
    //! Returns the grid origin.
    const vox::Vector3D& origin() const;
    
    //! Returns the grid spacing.
    const vox::Vector3D& gridSpacing() const;
    
    //! Returns the bounding box of the grid.
    const vox::BoundingBox3D& boundingBox() const;
    
    //! Returns the function that maps grid index to the cell-center position.
    DataPositionFunc cellCenterPosition() const;
    
    //!
    //! \brief Invokes the given function \p func for each grid cell.
    //!
    //! This function invokes the given function object \p func for each grid
    //! cell in serial manner. The input parameters are i, j, and k indices of a
    //! grid cell. The order of execution is i-first, j-next, k-last.
    //!
    void forEachCellIndex(
                          const std::function<void(uint, uint, uint)>& func) const;
    
    //!
    //! \brief Invokes the given function \p func for each grid cell parallelly.
    //!
    //! This function invokes the given function object \p func for each grid
    //! cell in parallel manner. The input parameters are i, j, and k indices of
    //! a grid cell. The order of execution can be arbitrary since it's
    //! multi-threaded.
    //!
    void parallelForEachCellIndex(
                                  const std::function<void(uint, uint, uint)>& func) const;
    
    //! Returns true if resolution, grid-spacing and origin are same.
    bool hasSameShape(const Grid3& other) const;
    
    //! Swaps the data with other grid.
    virtual void swap(Grid3* other) = 0;
    
protected:
    //! Sets the size parameters including the resolution, grid spacing, and
    //! origin.
    void setSizeParameters(const vox::Size3& resolution,
                           const vox::Vector3D& gridSpacing,
                           const vox::Vector3D& origin);
    
    //! Swaps the size parameters with given grid \p other.
    void swapGrid(Grid3* other);
    
    //! Sets the size parameters with given grid \p other.
    void setGrid(const Grid3& other);
    
private:
    vox::Size3 _resolution;
    vox::Vector3D _gridSpacing = vox::Vector3D(1, 1, 1);
    vox::Vector3D _origin;
    vox::BoundingBox3D _boundingBox = vox::BoundingBox3D(vox::Vector3D(),
                                                         vox::Vector3D());
};

typedef std::shared_ptr<Grid3> Grid3Ptr;

#define JET_GRID3_TYPE_NAME(DerivedClassName) \
std::string typeName() const override { return #DerivedClassName; }

}  // namespace vox

#endif  // INCLUDE_VDB_GRID3_H_
