////
////  openvdb_test.cpp
////  unittest
////
////  Created by Feng Yang on 2020/1/2.
////  Copyright © 2020 Feng Yang. All rights reserved.
////
//
//#include <openvdb/io/Stream.h>
//#include <openvdb/openvdb.h>
//
//#include "../external/gtest/include/gtest/gtest.h"
//
//TEST(GridGenerator3, createCellCenteredScalarGrid) {
//    openvdb::DoubleGrid::Ptr grid = openvdb::DoubleGrid::create();
//    // Get a voxel accessor.
//    typename openvdb::DoubleGrid::Accessor accessor = grid->getAccessor();
//    openvdb::Coord ijk;
//    int &i = ijk[0], &j = ijk[1], &k = ijk[2];
//    for (i = 0; i < 5; ++i) {
//        for (j = 0; j < 5; ++j) {
//            for (k = 0; k < 5; ++k) {
//                accessor.setValue(ijk, 0.0);
//            }
//        }
//    }
//
//    for (openvdb::DoubleGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter) {
//        std::cout<<iter.getCoord()<<"\n";
//    }
//    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//
//    using IterRange = openvdb::tree::IteratorRange<openvdb::DoubleGrid::ValueOnIter>;
//    IterRange range(grid->beginValueOn());
//    tbb::parallel_for(range, [&](IterRange& r) {
//        // Iterate over a subrange of the leaf iterator's iteration space.
//        for ( ; r; ++r) {
//            openvdb::DoubleGrid::ValueOnIter iter = r.iterator();
//            std::cout<<iter.getCoord()<<"\n";
//        }
//    });
//    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
//
//    using IterRange2 = openvdb::tree::IteratorRange<openvdb::DoubleGrid::TreeType::LeafIter>;
//    IterRange2 range2(grid->tree().beginLeaf());
//    tbb::parallel_for(range2, [&](IterRange2& r) {
//        // Iterate over a subrange of the leaf iterator's iteration space.
//        auto iter = r.iterator()->beginValueOn();
//        for ( ; iter; ++iter) {
//            std::cout<<iter.getCoord()<<"\n";
//        }
//    });
//}
