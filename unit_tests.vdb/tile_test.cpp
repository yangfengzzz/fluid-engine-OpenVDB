//
//  tile_test.cpp
//  vdb_tests
//
//  Created by Feng Yang on 2020/3/30.
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

/*
TEST(Tile, Level1) {
    openvdb::DoubleGrid::Ptr grid = openvdb::DoubleGrid::create(0.0);
    grid->treePtr()->addTile(1, openvdb::Coord(0, 0, 0), 10, true);
    auto iter = grid->treePtr()->beginValueOn();
    for (; iter; ++iter) {
        if (iter.isTileValue() == true) {
            std::cout<<"TRUE"<<"\n";
        }else{
            std::cout<<"FALSE"<<"\n";
        }
        std::cout<<iter.getCoord()<<"\n";
    }
    
    std::cout<<grid->treePtr()->getValue(openvdb::Coord(1, 0, 0))<<"------\n";
    grid->treePtr()->setValueOn(openvdb::Coord(1, 0, 0), 23);
    iter = grid->treePtr()->beginValueOn();
    for (; iter; ++iter) {
        if (iter.isTileValue() == true) {
            std::cout<<"TRUE"<<"\n";
        }else{
            std::cout<<"FALSE"<<"\n";
        }
        std::cout<<iter.getCoord()<<"\t"<<iter.getValue()<<"\n";
    }
    
    std::cout<<"======================\n";
    grid->treePtr()->addTile(1, openvdb::Coord(0, 0, 0), 20, true);
    iter = grid->treePtr()->beginValueOn();
    for (; iter; ++iter) {
        if (iter.isTileValue() == true) {
            std::cout<<"TRUE"<<"\n";
        }else{
            std::cout<<"FALSE"<<"\n";
        }
        std::cout<<iter.getCoord()<<"\t"<<iter.getValue()<<"\n";
    }
    std::cout<<grid->treePtr()->getValue(openvdb::Coord(1, 0, 0))<<"------\n";
}

TEST(Tile, Level2) {
    openvdb::DoubleGrid::Ptr grid = openvdb::DoubleGrid::create(0.0);
    grid->treePtr()->addTile(1, openvdb::Coord(0, 0, 0), 10, true);
    auto iter = grid->treePtr()->beginValueOn();
    for (; iter; ++iter) {
        if (iter.isTileValue() == true) {
            std::cout<<"TRUE"<<"\n";
        }else{
            std::cout<<"FALSE"<<"\n";
        }
        std::cout<<iter.getCoord()<<"\n";
    }
    
    grid->treePtr()->addTile(2, openvdb::Coord(0, 0, 0), 20, true);
    iter = grid->treePtr()->beginValueOn();
    for (; iter; ++iter) {
        if (iter.isTileValue() == true) {
            std::cout<<"TRUE"<<"\n";
        }else{
            std::cout<<"FALSE"<<"\n";
        }
        std::cout<<iter.getCoord()<<"\n";
    }
}
*/

TEST(Tile, Adaptive) {
    openvdb::DoubleGrid::Ptr grid = openvdb::DoubleGrid::create(0.0);
    grid->treePtr()->setValueOn(openvdb::Coord(0, 0, 0), 10);
    grid->treePtr()->addTile(1, openvdb::Coord(8, 0, 0), 20, true);
    auto iter = grid->treePtr()->beginValueOn();
    for (; iter; ++iter) {
        if (iter.isTileValue() == true) {
            std::cout<<"TRUE"<<"\n";
        }else{
            std::cout<<"FALSE"<<"\n";
        }
        std::cout<<iter.getCoord()<<"\n";
    }
    std::cout<<"~~~~~~~~~~\n";
    
    typedef openvdb::DoubleTree::RootNodeType RootT;
    typedef RootT::ChildNodeType Int1T;
    typedef Int1T::ChildNodeType Int2T;
    typedef Int2T::ChildNodeType LeafT;
    RootT* root; Int1T* int1; Int2T* int2; LeafT* leaf;
    
    iter = grid->treePtr()->beginValueOn();
    for (; iter; ++iter) {
        const openvdb::Index iterLevel = iter.getLevel();
        std::cout<<"LEVEL: "<<iterLevel<<"\n";
        iter.getNode(root);
        if (root != nullptr) {
            auto node_iter = root->beginValueOn();
            for (; node_iter; ++node_iter) {
                std::cout<<node_iter.getCoord()<<"---\n";
            }
        }

        iter.getNode(int1);
        if (int1 != nullptr) {
            auto node_iter = int1->beginValueOn();
            for (; node_iter; ++node_iter) {
                std::cout<<node_iter.getCoord()<<"===\n";
            }
        }
        
        iter.getNode(int2);
        if (int2 != nullptr) {
            auto node_iter = int2->beginValueOn();
            for (; node_iter; ++node_iter) {
                std::cout<<node_iter.getCoord()<<"~~~\n";
            }
        }
        
        iter.getNode(leaf);
        if (leaf != nullptr) {
            auto node_iter = leaf->beginValueOn();
            for (; node_iter; ++node_iter) {
                std::cout<<node_iter.getCoord()<<"@@@\n";
            }
        }
        std::cout<<"++++++++++++++++++++\n";
    }
    
}
