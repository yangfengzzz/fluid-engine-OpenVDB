//
//  point_VDB_searcher3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/1/5.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef point_VDB_searcher3_hpp
#define point_VDB_searcher3_hpp

#include "../src.common/point_neighbor_searcher3.h"

#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <vector>

namespace vdb {

class PointList
{
public:
    typedef openvdb::Vec3R  PosType;
    
    PointList(const std::vector<PosType>& points)
    : mPoints(&points)
    {
    }
    
    size_t size() const {
        return mPoints->size();
    }
    
    void getPos(size_t n, PosType& xyz) const {
        xyz = (*mPoints)[n];
    }
    
protected:
    std::vector<PosType> const * const mPoints;
}; // PointList

class MyParticleList
{
protected:
    struct MyParticle {
        openvdb::Vec3R p, v;
        openvdb::Real  r;
    };
    openvdb::Real           mRadiusScale;
    openvdb::Real           mVelocityScale;
    std::vector<MyParticle> mParticleList;
public:
    
    typedef openvdb::Vec3R  PosType;
    
    MyParticleList(openvdb::Real rScale=1, openvdb::Real vScale=1)
    : mRadiusScale(rScale), mVelocityScale(vScale) {}
    void add(const openvdb::Vec3R &p, const openvdb::Real &r,
             const openvdb::Vec3R &v=openvdb::Vec3R(0,0,0))
    {
        MyParticle pa;
        pa.p = p;
        pa.r = r;
        pa.v = v;
        mParticleList.push_back(pa);
    }
    /// @return coordinate bbox in the space of the specified transfrom
    openvdb::CoordBBox getBBox(const openvdb::GridBase& grid) {
        openvdb::CoordBBox bbox;
        openvdb::Coord &min= bbox.min(), &max = bbox.max();
        openvdb::Vec3R pos;
        openvdb::Real rad, invDx = 1/grid.voxelSize()[0];
        for (size_t n=0, e=this->size(); n<e; ++n) {
            this->getPosRad(n, pos, rad);
            const openvdb::Vec3d xyz = grid.worldToIndex(pos);
            const openvdb::Real   r  = rad * invDx;
            for (int i=0; i<3; ++i) {
                min[i] = openvdb::math::Min(min[i], openvdb::math::Floor(xyz[i] - r));
                max[i] = openvdb::math::Max(max[i], openvdb::math::Ceil( xyz[i] + r));
            }
        }
        return bbox;
    }
    //typedef int AttributeType;
    // The methods below are only required for the unit-tests
    openvdb::Vec3R pos(int n)   const {return mParticleList[n].p;}
    openvdb::Vec3R vel(int n)   const {return mVelocityScale*mParticleList[n].v;}
    openvdb::Real radius(int n) const {return mRadiusScale*mParticleList[n].r;}
    
    //////////////////////////////////////////////////////////////////////////////
    /// The methods below are the only ones required by tools::ParticleToLevelSet
    /// @note We return by value since the radius and velocities are modified
    /// by the scaling factors! Also these methods are all assumed to
    /// be thread-safe.
    
    /// Return the total number of particles in list.
    ///  Always required!
    size_t size() const { return mParticleList.size(); }
    
    /// Get the world space position of n'th particle.
    /// Required by ParticledToLevelSet::rasterizeSphere(*this,radius).
    void getPos(size_t n,  openvdb::Vec3R&pos) const { pos = mParticleList[n].p; }
    
    
    void getPosRad(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
    }
    void getPosRadVel(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad, openvdb::Vec3R& vel) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
        vel = mVelocityScale*mParticleList[n].v;
    }
    // The method below is only required for attribute transfer
    void getAtt(size_t n, openvdb::Index32& att) const { att = openvdb::Index32(n); }
};

//!
//! \brief VDB 3-D point searcher.
//!
//! This class implements 3-D point searcher simply by looking up every point in
//! the VDB.
//!
class PointVDBSearcher3 final : public vox::PointNeighborSearcher3 {
public:
    JET_NEIGHBOR_SEARCHER3_TYPE_NAME(PointVDBSearcher3)
    
    class Builder;
    
    //! Default constructor.
    PointVDBSearcher3();
    
    //! Copy constructor.
    PointVDBSearcher3(const PointVDBSearcher3& other);
    
public:
    //! Builds internal acceleration structure for given points list.
    void build(const vox::ConstArrayAccessor1<vox::Vector3D>& points) override;
    
    //!
    //! Invokes the callback function for each nearby point around the origin
    //! within given radius.
    //!
    //! \param[in]  origin   The origin position.
    //! \param[in]  radius   The search radius.
    //! \param[in]  callback The callback function.
    //!
    void forEachNearbyPoint(
                            const vox::Vector3D& origin, double radius,
                            const ForEachNearbyPointFunc& callback) override;
    
    //!
    //! Returns true if there are any nearby points for given origin within
    //! radius.
    //!
    //! \param[in]  origin The origin.
    //! \param[in]  radius The radius.
    //!
    //! \return     True if has nearby point, false otherwise.
    //!
    bool hasNearbyPoint(const vox::Vector3D& origin, double radius) override;
    
    openvdb::tools::PointIndexGrid::Ptr getIndexGrid();
    
    openvdb::points::PointDataGrid::Ptr getDataGrid();
    
    openvdb::math::Transform::Ptr getTransform();
    
    openvdb::Vec3d getPoint(size_t i);
    
public:
    //!
    //! \brief      Creates a new instance of the object with same properties
    //!             than original.
    //!
    //! \return     Copy of this object.
    //!
    vox::PointNeighborSearcher3Ptr clone() const override;
    
    //! Assignment operator.
    PointVDBSearcher3& operator=(const PointVDBSearcher3& other);
    
    //! Copy from the other instance.
    void set(const PointVDBSearcher3& other);
    
    //! Serializes the neighbor searcher into the buffer.
    void serialize(std::vector<uint8_t>* buffer) const override;
    
    //! Deserializes the neighbor searcher from the buffer.
    void deserialize(const std::vector<uint8_t>& buffer) override;
    
    //! Returns builder fox PointVDBSearcher3.
    static Builder builder();
    
private:
    std::vector<openvdb::Vec3R> _points;
    
    openvdb::math::Transform::Ptr transform;
    
    openvdb::tools::PointIndexGrid::Ptr pointIndexGrid;
    
    openvdb::points::PointDataGrid::Ptr pointDataGrid;
    
    openvdb::tools::PointIndexIterator<> iterator;
};

//! Shared pointer for the PointSimpleVDBSearcher3 type.
typedef std::shared_ptr<PointVDBSearcher3> PointVDBSearcher3Ptr;

//!
//! \brief Front-end to create PointSimpleVDBSearcher3 objects step by step.
//!
class PointVDBSearcher3::Builder final
: public vox::PointNeighborSearcherBuilder3 {
public:
    //! Builds PointVDBSearcher3 instance.
    PointVDBSearcher3 build() const;
    
    //! Builds shared pointer of PointVDBSearcher3 instance.
    PointVDBSearcher3Ptr makeShared() const;
    
    //! Returns shared pointer of PointNeighborSearcher3 type.
    vox::PointNeighborSearcher3Ptr buildPointNeighborSearcher() const override;
};

}  // namespace vox

#endif /* point_VDB_searcher3_hpp */
