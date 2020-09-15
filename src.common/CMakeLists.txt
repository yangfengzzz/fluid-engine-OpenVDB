cmake_minimum_required(VERSION 3.14)

aux_source_directory(./ DIR_SRCS)

list(APPEND CMAKE_MODULE_PATH "../cmake")
find_package(glog)
find_package(TBB)

set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -fPIC -g ${CMAKE_CXX_FLAGS}")
set(CMAKE_LD_FLAGS "-ldl ${CMAKE_LD_FLAGS}")

ADD_LIBRARY(Flex_Common SHARED ${DIR_SRCS})   
target_link_libraries(Flex_Common glog::glog TBB::tbb)
