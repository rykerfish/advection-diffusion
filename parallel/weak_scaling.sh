#!/bin/bash
#
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --time=01:00:00

echo "Running weak scaling test"
echo "Job ran on: $SLURM_NODELIST"

n_base=16

module load mpi/openmpi/gcc/4.1.1
mpicc -o main main.c parallel.c -lm

for((i=0; i<=2; i+=1)); do

    p=$((4**i))
    n=$((n_base*2**i))
    echo "Testing with $p threads."

    time mpirun -np $p ./main $n

done

echo "Job finished."
