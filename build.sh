#!/bin/bash

# Move to build directory
cd code_students/build

# load module for mpi
module load gcc/12.2.0-gcc-8.5.0-p4pe45v openmpi/3.1.6-gcc-12.2.0-d2gmn55

# Generate Makefiles using CMake, using hdf in phillips scratch dir
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="/scratch/c703429/software/hdf5-1.14.4-2" ..

# Build project using Make
make
