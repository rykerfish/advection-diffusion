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

	// declare vectors
	Vector u_grid; 
	Vector x_grid;
	Vector y_grid;
	
	// allocate vectors
	int rows = grid_size;
	int cols = grid_size;
	allocate_vector(&u_grid, rows, cols);
	allocate_vector(&x_grid, rows, cols);
	allocate_vector(&y_grid, rows, cols);
	
	// create grids for constructing the concentration grid
	float x_start = -3.0;
	float x_end = 3.0;
	float y_start = -3.0;
	float y_end = 3.0;
	create_grid(&x_grid, x_start, x_end, 'x');
	create_grid(&y_grid, y_start, y_end, 'y');

	float dx = x_grid.data[1] - x_grid.data[0];
	float dy = y_grid.data[rows + 1] - y_grid.data[0];

	// construct the concentration grid
	initialize_concentration_vector(&u_grid, &x_grid, &y_grid);
	
	// determine current coordinates
	int rank;
	int coords[ndims];
	MPI_Comm_rank(grid, &rank);
	MPI_Cart_coords(grid, rank, ndims, coords);
	
	// calculate size of subgrid
	sub_x_dim = u_grid.x_dim / dims[0];
	sub_y_dim = u_grid.y_dim / dims[1];
	
	// determine start and end indices
	int x_start = coords[0] * sub_x_dim;
	int y_start = coords[1] * sub_y_dim;
	int x_end = x_start + sub_x_dim;
	int y_end = y_start + sub_y_dim;
	
	// make sure these are valid indices
	if(x_end > u_grid.x_dim || y_end > u_grid.y_dim){
		printf("Error: sub arrays are not contained within main array");
		return 1;
	}
	
	// allocated memory for local array
	float* local data = (float*)malloc(sub_x_dim * sub_y_dim * sizeof(float));
	
	// copy the relevant data
	int i, j;
	for(i = x_start; i < x_end; ++i) {
		for(j = y_start; j < y_end; ++j){
			int sub_i = i - x_start;
			int sub_j = j - y_start;
			int idx = i * u_grid.y_dim + j;
			int sub_idx = sub_i * sub_y_dim + sub_j;
			local_data[sub_idx] = u_grid.data[idx];
		}
	}
	

}
