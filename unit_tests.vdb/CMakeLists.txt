CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

project(unit_tests)
aux_source_directory(./ DIR_SRCS)

list(APPEND CMAKE_MODULE_PATH "/Users/yangfeng/Desktop/DigitalFlex2/cmake")
find_package(OpenVDB REQUIRED)
find_package(TBB)
find_package(glog)

include_directories("../external/gtest/include")
include_directories("../external/gmock/include")
link_directories("../external/gtest")
link_directories("../external/gmock")
link_directories("../build/src.common")
link_directories("../build/src.vdb")

set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -fPIC -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_LD_FLAGS "${CMAKE_LD_FLAGS}")

add_executable(app ${DIR_SRCS})
target_link_libraries(app gtest_main gmock_main glog::glog TBB::tbb OpenVDB::openvdb Flex_Common Flex_VDB)
