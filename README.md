Please Refer to [DigitalArche2](https://github.com/yangfengzzz/DigitalArche2) for the further development which will include rendering, Skeleton and GPU simulations.
----------------------------------------------
This repos is inspired by [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev) which teach me how to code numerical algorithm by using C++. In this repos, files in src.common folder are all forked from [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev). And all my own work is located in src.vdb, unit_tests.vdb and manual_tests.vdb, which rewrite 3D Grids by using OpenVDB.

[openvdb](https://github.com/AcademySoftwareFoundation/openvdb) is a well-known opensource C++ code which implemented the sparse data structure called VDB. In  the physical simualtion of computer graphics, sparsity can really save the memory and increase the speed of simualtions. This is the main reason I want to implement this repos when I learned from [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev).

My own code is not well optimized by using all features of the project OpenVDB, which not only contain VDB but a bunch of algorithm about levelset treatment. This project is just a attempt which can help me to better understant [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev). But I also think this repos can inspire some future works which use more functionality of sparse data structure.

# INSTALL

This repos is depended on all dependence of [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev) including TBB and OpenVDB. You can refer [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev) to install the dependence which is very easy in Mac by using Homebrew. 

After installing all the dependences, you can just run the script:

```bash
./install
```

which will call cmake and compile the source code by using 4 core (which you can adjust in the install.sh)

# USAGE

You can choose hybrid method like FLIP in manualTests.vdb and import the results into vdb files. Houdini support VDB and render the sparse grid very well. you can find the video in my [bilibili](https://www.bilibili.com/video/BV1XV41127Cy/) channel which will show the Dam Break Problem by using FLIP.
