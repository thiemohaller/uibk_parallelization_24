#!/bin/bash

# Help menu
case $1 in
    -h|--help)
        echo "Usage: $0 [OPTIONS]"
        echo "Options:"
        echo "  -h, --help     Show help menu"
        echo "  -b, --build    Build the project before submitting"
        exit 0
        ;;
    -b|--build)
        echo "🔨 Building the project"
        ./build.sh
        ;;
esac

echo "🚀 Submitting parallel job via slurm"
sbatch job_parallel.sh

echo "👀 Watching submitted job"
watch -n 1 squeue --me

# echo "📂 Script location: /scratch/cb761148/uibk_parallelization_24/run_parallel.sh"
