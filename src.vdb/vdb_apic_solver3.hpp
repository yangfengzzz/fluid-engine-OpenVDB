//
//  vdb_apic_solver3.hpp
//  Solvers
//
//  Created by Feng Yang on 2020/2/13.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#ifndef INCLUDE_VDB_APIC_SOLVER3_H_
#define INCLUDE_VDB_APIC_SOLVER3_H_

#include "vdb_pic_solver3.hpp"

namespace vdb {

//!
//! \brief 3-D Affine Particle-in-Cell (APIC) implementation
//!
//! This class implements 3-D Affine Particle-in-Cell (APIC) solver from the
//! SIGGRAPH paper, Jiang 2015.
//!
//! \see Jiang, Chenfanfu, et al. "The affine particle-in-cell method."
//!      ACM Transactions on Graphics (TOG) 34.4 (2015): 51.
//!
class ApicSolver3 : public PicSolver3 {
public:
    class Builder;
    
    //! Default constructor.
    ApicSolver3();
    
    //! Constructs solver with initial grid size.
    ApicSolver3(const vox::Size3& resolution,
                const vox::Vector3D& gridSpacing,
                const vox::Vector3D& gridOrigin);
    
    //! Default destructor.
    virtual ~ApicSolver3();
    
    //! Returns builder fox ApicSolver3.
    static Builder builder();
    
protected:
    //! Transfers velocity field from particles to grids.
    void transferFromParticlesToGrids() override;
    
    //! Transfers velocity field from grids to particles.
    void transferFromGridsToParticles() override;
    
private:
    vox::Array1<vox::Vector3D> _cX;
    vox::Array1<vox::Vector3D> _cY;
    vox::Array1<vox::Vector3D> _cZ;
};

//! Shared pointer type for the ApicSolver3.
typedef std::shared_ptr<ApicSolver3> ApicSolver3Ptr;


//!
//! \brief Front-end to create ApicSolver3 objects step by step.
//!
class ApicSolver3::Builder final
: public GridFluidSolverBuilderBase3<ApicSolver3::Builder> {
public:
    //! Builds ApicSolver3.
    ApicSolver3 build() const;
    
    //! Builds shared pointer of ApicSolver3 instance.
    ApicSolver3Ptr makeShared() const;
};

}  // namespace vox

#endif  // INCLUDE_VDB_APIC_SOLVER3_H_

