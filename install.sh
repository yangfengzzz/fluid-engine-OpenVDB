#!/bin/bash

#compile the shared library
echo "compile the shared library--------------------------"
mkdir -pv build
cd build
cmake ..
make -j4
echo "finish ----------------------------------------------"
