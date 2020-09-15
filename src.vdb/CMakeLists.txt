cmake_minimum_required(VERSION 3.14)

aux_source_directory(./ DIR_SRCS)

list(APPEND CMAKE_MODULE_PATH "../cmake")
find_package(OpenVDB REQUIRED)
find_package(TBB)
find_package(glog)

set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -fPIC -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_LD_FLAGS "-ldl ${CMAKE_LD_FLAGS}")

ADD_LIBRARY(Flex_VDB SHARED ${DIR_SRCS})   
target_link_libraries(Flex_VDB Flex_Common OpenVDB::openvdb glog::glog TBB::tbb)
