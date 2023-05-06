buildMain:
	mpicc -g -o main main.c parallel.c -lm 

runMain:
	./main

cleanMain:
	rm main

test:
	gcc -g -o test test.c parallel.c -lm
	./test
	rm ./test

test_par:
	module load mpi/openmpi/gcc/4.1.1
	mpicc -g -o test_par test_parallel.c parallel.c -lm
	mpirun -np 3 ./test_par
	rm ./test_par