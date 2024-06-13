#!/bin/bash

echo "Submitting mandelbrot job via slurm"
sbatch job.sh

echo "Watching submitted job"
# watch -n 1 squ
watch -n 1 squeue --me
# squ is an alias for squeue --me
