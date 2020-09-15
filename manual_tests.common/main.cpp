// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "manual_tests.h"
#include "../external/gtest/include/gtest/gtest.h"
#include <fstream>

using namespace vox;

int main(int argc, char** argv) {    
    ::testing::InitGoogleTest(&argc, argv);
    createDirectory(JET_TESTS_OUTPUT_DIR);
    
    int ret = RUN_ALL_TESTS();
    
    return ret;
}
