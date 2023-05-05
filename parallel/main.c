#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parallel.h"

int main(int argc, char** argv) {

    if(argc < 2){
        printf("Error: Should call with <# of processes> <grid size>\n");
        return 1;
    }
    int n_procs = atoi(argv[1]);
    int grid_size = atoi(argv[2]);

    MPI_Init(&argc, &argv);

    int ndims = 2;
    int dims[ndims];
	int periodic[] = {1, 1};
    MPI_Dims_create(n_procs, ndims, dims);

	MPI_Comm grid;	
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, 1, &grid);

	// // declare vectors
	// Vector u_grid; 
	// Vector x_grid;
	// Vector y_grid;
	
	// // allocate vectors
	// int rows = grid_size;
	// int cols = grid_size;
	// allocate_vector(&u_grid, rows, cols);
	// allocate_vector(&x_grid, rows, cols);
	// allocate_vector(&y_grid, rows, cols);
	
	// // create grids for constructing the concentration grid
	// float x_start = -3.0;
	// float x_end = 3.0;
	// float y_start = -3.0;
	// float y_end = 3.0;
	// create_grid(&x_grid, x_start, x_end, 'x');
	// create_grid(&y_grid, y_start, y_end, 'y');

	// float dx = x_grid.data[1] - x_grid.data[0];
	// float dy = y_grid.data[rows + 1] - y_grid.data[0];

	// // construct the concentration grid
	// initialize_concentration_vector(&u_grid, &x_grid, &y_grid);





}