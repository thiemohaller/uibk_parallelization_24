#!/bin/bash

BUILD=false
DEBUG=false

# Parse arguments to get all added flags
while [[ "$1" != "" ]]; do
    case $1 in
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  -h, --help     Show help menu"
            echo "  -b, --build    Build the project before submitting"
            echo "  -d, --debug    Submit debug job instead of regular job"
            exit 0
            ;;
        -b|--build)
            BUILD=true
            ;;
        -d|--debug)
            DEBUG=true
            ;;
        *)
            echo "Invalid option: $1"
            exit 1
            ;;
    esac
    shift
done

# rebuild project (calls cmake and make)
if [ "$BUILD" = true ]; then
    echo "🔨 Building the project"
    ./build.sh
fi

# submit debug version using gdb
if [ "$DEBUG" = true ]; then
    echo "🤖 Submitting debug job via slurm, using valgrind to trace leaks, excellent choice"
    sbatch job_parallel_debug.sh
else
    echo "🚀 Submitting regular job via slurm"
    sbatch job_parallel.sh
fi

echo "👀 Use \`watch -n 1 squeue --me\` to watch your submitted jobs"
echo "🐬 So long, and thanks for all the fish! 🐬"
