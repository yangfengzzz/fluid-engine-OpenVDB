//
//  IndexTree_tests.cpp
//  vdb_tests
//
//  Created by Feng Yang on 2020/2/11.
//  Copyright © 2020 Feng Yang. All rights reserved.
//

#include "../src.vdb/vdb_cell_centered_scalar_grid3.h"
#include "../src.vdb/vdb_face_centered_grid3.h"
#include "../src.vdb/vdb_helper.h"
#include <openvdb/tools/PoissonSolver.h>
#include <vector>
#include <iterator>

#include "../external/gtest/include/gtest/gtest.h"

using namespace vdb;

TEST(IndexTree, Construction) {
    openvdb::DoubleGrid::Ptr finest = openvdb::DoubleGrid::create(0.0);
    openvdb::math::ScaleTranslateMap transMap(openvdb::Vec3d(0.001, 0.001, 0.001),
                                              openvdb::Vec3d(1, 2, 3));
    openvdb::math::Transform trans(transMap.copy());
    finest->setTransform(trans.copy());
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                finest->tree().setValueOn(openvdb::Coord(i, j, k), 100);
            }
        }
    }
    
    using VIdxTreeT = openvdb::DoubleTree::template ValueConverter<openvdb::tools::poisson::VIndex>::Type;
    typename VIdxTreeT::Ptr idxTree = openvdb::tools::poisson::createIndexTree(finest->tree());
    
    std::vector<int> idx;
    for (VIdxTreeT::ValueOnIter iter = idxTree->beginValueOn(); iter; ++iter) {
        idx.push_back(idxTree->getValue(iter.getCoord()) );
    }
    EXPECT_EQ(idx.size(), 1000);
    std::sort(idx.begin(), idx.end());
    for (int i = 0; i < idx.size(); i++) {
        EXPECT_EQ(idx[i], i);
    }
}

TEST(IndexTree, createVectorFromTree) {
    openvdb::DoubleGrid::Ptr finest = openvdb::DoubleGrid::create(0.0);
    openvdb::math::ScaleTranslateMap transMap(openvdb::Vec3d(0.001, 0.001, 0.001),
                                              openvdb::Vec3d(1, 2, 3));
    openvdb::math::Transform trans(transMap.copy());
    finest->setTransform(trans.copy());
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                finest->tree().setValueOn(openvdb::Coord(i, j, k), i+j+k);
            }
        }
    }
    
    using VIdxTreeT = openvdb::DoubleTree::template ValueConverter<openvdb::tools::poisson::VIndex>::Type;
    typename VIdxTreeT::ConstPtr
    idxTree = openvdb::tools::poisson::createIndexTree(finest->tree());
    
    openvdb::math::pcg::Vector<double>::Ptr
    b = openvdb::tools::poisson::createVectorFromTree<double>(finest->tree(),
                                                              *idxTree);
    
    std::vector<double> rhs(b->data(), b->data()+b->size());
    //    for (int i = 0; i < rhs.size(); i++) {
    //        std::cout<<rhs[i]<<"\t"<<b->at(i)<<"\n";
    //    }
}

TEST(IndexTree, createTreeFromVector) {
    openvdb::DoubleGrid::Ptr finest = openvdb::DoubleGrid::create(0.0);
    openvdb::math::ScaleTranslateMap transMap(openvdb::Vec3d(0.001, 0.001, 0.001),
                                              openvdb::Vec3d(1, 2, 3));
    openvdb::math::Transform trans(transMap.copy());
    finest->setTransform(trans.copy());
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                finest->tree().setValueOn(openvdb::Coord(i, j, k), i+j+k);
            }
        }
    }
    
    using VIdxTreeT = openvdb::DoubleTree::template ValueConverter<openvdb::tools::poisson::VIndex>::Type;
    typename VIdxTreeT::ConstPtr
    idxTree = openvdb::tools::poisson::createIndexTree(finest->tree());
    
    openvdb::math::pcg::Vector<double>::Ptr
    b = openvdb::tools::poisson::createVectorFromTree<double>(finest->tree(),
                                                              *idxTree);
    
    openvdb::DoubleTree::Ptr
    back = openvdb::tools::poisson::createTreeFromVector<double>(*b, *idxTree, 0.0);
    for (openvdb::DoubleGrid::ValueOnIter iter = finest->beginValueOn(); iter; ++iter) {
        EXPECT_EQ(finest->tree().getValue(iter.getCoord()), back->getValue(iter.getCoord()));
    }
}
