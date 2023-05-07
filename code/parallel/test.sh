#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:5:00


#Record the node that we ran on
echo "RUNNING TESTS ON: $SLURM_NODELIST"

#Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# run test suite
mpicc -o test_par test_parallel.c parallel.c -lm
mpirun -np 3 test_par

# run full sim test
mpicc -o main main.c parallel.c -lm
mpirun -np 4 main 128
echo "Wrote output file from full test run to test_out/parallel_sim_out.txt"

echo "PARALLEL TESTS FINISHED"
