This repos is inspired by [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev) which teach me how to code numerical algorithm by using C++. In this repos, files in src.common folder are all forked from [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev). And all my own work is located in src.vdb, unit_tests.vdb and manual_tests.vdb, which rewrite 3D Grids by using OpenVDB.

[openvdb](https://github.com/AcademySoftwareFoundation/openvdb) is a well-known opensource C++ code which implemented the sparse data structure called VDB. In  the physical simualtion of computer graphics, sparsity can really save the memory and increase the speed of simualtions. This is the main reason I want to implement this repos when I learned from [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev).

My own code is not well optimized by using all features of the project OpenVDB, which not only contain VDB but a bunch of algorithm about levelset treatment. This project is just a attempt which can help me to better understant [fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev). But I also think this repos can inspire some future works which use more functionality of sparse data structure.