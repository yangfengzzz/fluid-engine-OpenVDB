#!/bin/bash

mkdir build
cd build
cmake -Wno-dev ..
make -j6
cd ..
ln -sf build/app app
