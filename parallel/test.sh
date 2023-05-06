#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=3
#SBATCH -t 00:01:00


#Record the node that we ran on
echo "RUNNING TESTS ON: $SLURM_NODELIST"

#Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

#Build code
mpicc -o test_par test_parallel.c parallel.c -lm

mpirun -np 3 test_par
