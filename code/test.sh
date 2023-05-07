#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:10:00

cd ./serial
echo "RUNNING SERIAL TESTS"
sbatch --output=../slurm_out/serial_test_results.txt ./test.sh
cd ..

cd ./parallel
echo "RUNNING PARALLEL TESTS"
sbatch --output=../slurm_out/parallel_test_results.txt ./test.sh
cd ..