#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:05:00

echo "RUNNING TESTS ON: $SLURM_NODELIST"

gcc -g -o test test.c serial.c -lm
./test
rm ./test

gcc -g -o main main.c serial.c -lm 
./main 128
echo "Wrote output file from full test run to test_out/serial_sim_out.txt"
rm ./main

echo "SERIAL TESTS FINISHED"