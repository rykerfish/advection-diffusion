#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:10:00

#Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

#Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

#Build code
mpicc -o main main.c parallel.c -lm

mpirun -np 4 main 128

echo "JOB FINISHED"