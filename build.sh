#!/bin/bash

# Default directory
directory="code_parallel"

# Function to display help menu
display_help() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -s, --sequential    Use code_students directory"
    echo "  -p, --parallel      Use code_parallel directory (default)"
    echo "  -h, --help          Display this help menu"
    exit 0
}

# Parse command-line options
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sequential)
            directory="code_students"
            shift
            ;;
        -p|--parallel)
            directory="code_parallel"
            shift
            ;;
        -h|--help)
            display_help
            ;;
        *)
            echo "Invalid option: $1" >&2
            exit 1
            ;;
    esac
done

# Move to build directory
cd "$directory/build"

# load module for mpi
module load gcc/12.2.0-gcc-8.5.0-p4pe45v openmpi/3.1.6-gcc-12.2.0-d2gmn55

# Generate Makefiles using CMake, using hdf in phillips scratch dir
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="/scratch/c703429/software/hdf5-1.14.4-2" ..

# Build project using Make
make
