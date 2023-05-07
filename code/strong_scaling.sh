#!/bin/bash
#
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --time=01:00:00

echo "Running strong scaling test"
echo "Job ran on: $SLURM_NODELIST"

n=256
max_procs=16

module load mpi/openmpi/gcc/4.1.1
mpicc -o main main.c parallel.c -lm

for((p=1; p<=$max_procs; p*=2)); do

    echo "Testing with $p threads."

    time mpirun -np $p main $n

done

echo "Job finished."