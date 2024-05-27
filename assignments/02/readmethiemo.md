# How to use?

## How to build

First time building: `mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make`

If build relevant files have changed, call cmake again, otherwise only calling `make` IN THE BUILD directory is enough.

## How to run

In root directory, in order to run with slurm, use `sbatch job.sh`.
Check if running with `squeue -u cb761148` or `squeue --me`.
