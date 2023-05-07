default:
	mpicc -g -o main main.c parallel.c -lm 

run:
	mpirun -np 4 ./main 128

clean:
	rm main

test:
	# module load mpi/openmpi/gcc/4.1.1
	mpicc -g -o test_par test_parallel.c parallel.c -lm
	mpirun -np 3 ./test_par
	rm ./test_par