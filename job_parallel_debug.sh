#!/bin/bash

# Execute job in the partition "lva" unless you have special requirements.
#SBATCH --partition=lva
# Name your job to be able to identify it later
#SBATCH --job-name test
# Redirect output stream to this file
#SBATCH --output=output_parallel.log
# Maximum number of tasks (=processes) to start in total
#SBATCH --ntasks=8
# Maximum number of tasks (=processes) to start per node
# SBATCH --ntasks-per-node=8
# Enforce exclusive node allocation, do not share with other jobs
#SBATCH --exclusive

module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

# /usr/bin/time -v mpiexec -n $SLURM_NTASKS ./code_parallel/build/apps/run_full_code.parallel
mpiexec -n 1 valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes ./code_parallel/build/apps/run_full_code.parallel :  -n 7 ./code_parallel/build/apps/run_full_code.parallel 2> output_parallel_valgrind.txt
