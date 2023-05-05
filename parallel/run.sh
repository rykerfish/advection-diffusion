#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:01:00

#Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

#Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

#Build code
mpicc -o main main.c parallel.c -lm

./main 4 100