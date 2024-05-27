#!/bin/bash

# Execute job in the partition "lva" unless you have special requirements.
#SBATCH --partition=lva
# Name your job to be able to identify it later
#SBATCH --job-name mandelbrot_tha
# Redirect output stream to this file
#SBATCH --output=mandelbrot.log
# Maximum number of tasks (=processes) to start in total
#SBATCH --ntasks=72
# Maximum number of tasks (=processes) to start per node, LLC3 has 12 cpus per core without hyperthreading
#SBATCH --ntasks-per-node=12
# Enforce exclusive node allocation, do not share with other jobs
#SBATCH --exclusive

module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

mpiexec -n $SLURM_NTASKS ./build/mandelbrot_mpi 15360 8640
