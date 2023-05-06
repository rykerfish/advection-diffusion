#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parallel.h"

int main(int argc, char** argv) {

	int i;

	if(argc < 1){
		printf("Error: Should call with <# of processes> <grid size>\n");
		return 1;
	}
	
	int grid_size = atoi(argv[1]);

	MPI_Init(&argc, &argv);
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// determine current rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	Matrix u_grid; 
	if(rank == 0){
		// declare matrixs
		Matrix x_grid;
		Matrix y_grid;

		// allocate matrixs
		int rows = grid_size;
		int cols = grid_size;
		allocate_matrix(&u_grid, rows, cols);
		allocate_matrix(&x_grid, rows, cols);
		allocate_matrix(&y_grid, rows, cols);

		// create grids for constructing the concentration grid
		float a = -3.0; 
		float b = 3.0;
		float x_bounds[] = {a, b};
		float y_bounds[] = {a, b};
		create_grid(&x_grid, x_bounds[0], x_bounds[1], 'x');
		create_grid(&y_grid, y_bounds[0], y_bounds[1], 'y');

		float dx = x_grid.data[1] - x_grid.data[0];
		float dy = y_grid.data[rows + 1] - y_grid.data[0];

		// construct the concentration grid
		initialize_concentration_matrix(&u_grid, &x_grid, &y_grid);

		deallocate_matrix(&x_grid);	
		deallocate_matrix(&y_grid);
	}

	
	// calculate size of subgrid
	int row_len = grid_size;
	int num_rows_in_block = grid_size / nprocs;
	int local_size = row_len * num_rows_in_block; 
	int padded_size = local_size + 2*row_len; // add 2 for ghost region rows

	// allocated memory for local array
	float* local_data = (float*)malloc(local_size * sizeof(float));
	
	// scatter the data to each process
	// start row_len into local_data so that there is room for ghost region
	MPI_Scatter(u_grid.data, local_size, MPI_FLOAT, local_data, local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	float* padded_data = (float*)malloc(padded_size * sizeof(float));
	for(i = 0; i <row_len; i++){
		padded_data[i+row_len] = local_data[i];
	}

	int up_neighbor, down_neighbor;
	if(rank == 0){
		up_neighbor = 3;
	} else {
		up_neighbor = rank - 1;
	}

	if(rank == 3){
		down_neighbor = 0;
	} else {
		down_neighbor = rank + 1;
	}
	
	// setup for communcations
	MPI_Status status[4];
	MPI_Request request[4];
	for(i = 0; i < 4; i++) request[i] = MPI_REQUEST_NULL;

	// do ghost region comms
	MPI_Isend(padded_data + row_len, row_len, MPI_FLOAT, up_neighbor, 0, MPI_COMM_WORLD, &request[0]); // send first row of data to upper neighbor
	printf("rank %d made 1\n", rank);
	MPI_Isend(padded_data + local_size, row_len, MPI_FLOAT, down_neighbor, 0, MPI_COMM_WORLD, &request[1]); // send last row of data to lower neighbor
	printf("rank %d made 2\n", rank);
	MPI_Irecv(padded_data, row_len, MPI_FLOAT, up_neighbor, 0, MPI_COMM_WORLD, &request[2]); // recv data from upper neighbor
	printf("rank %d made 3\n", rank);
	MPI_Irecv(padded_data + padded_size - row_len, row_len, MPI_FLOAT, down_neighbor, 0, MPI_COMM_WORLD, &request[3]); // recv data from lower neighbor
	printf("rank %d made 4\n", rank);

	MPI_Waitall(4, request, status);
	printf("DONE WAITING\n");


	// // set time step parameters
	// float cfl = 0.01;
	// float dt = cfl * dx;
	// float t_final = 10;
	// int n_steps = t_final / dt;

	// // set physical parameters
	// float diffusion = 0.01;

	// // do the simulation
	// Matrix u_update;
	// allocate_matrix(&u_update, rows, cols);

	// loop through the timesteps
	// int n;
	// for(n = 0; n < n_steps; ++n){
	// 	float u_lap, u_adv;

	// 	int x_i, y_i;
	// 	for(x_i = 0; x_i < cols; ++x_i){
	// 		for(y_i = 0; y_i < rows; ++y_i){
	// 			int idx = x_i + cols*y_i;
	// 			float u_center = u_grid.data[idx];
				
	// 			// 2D finite laplacian with finite boundary conditions
	// 			float u_west, u_east, u_north, u_south;
				
	// 			if(x_i == 0){ // left edge
	// 				u_west = u_grid.data[idx + cols - 1];
	// 			} else {
	// 				u_west = u_grid.data[idx - 1];
	// 			}

	// 			if(x_i == cols - 1){ // right edge
	// 				u_east = u_grid.data[idx - cols + 1];
	// 			} else {
	// 				u_east = u_grid.data[idx + 1];
	// 			}

	// 			if(y_i == 0){ // top edge
	// 				u_north = u_grid.data[(rows - 1)*cols + x_i];
	// 			} else {
	// 				u_north = u_grid.data[idx - cols];
	// 			}

	// 			if(y_i == rows - 1){ // bottom edge
	// 				u_south = u_grid.data[x_i];
	// 			} else {
	// 				u_south = u_grid.data[idx + cols];
	// 			}

	// 			// calculate finite difference laplacian with a 5 point stencil
	// 			u_lap = (u_west + u_east + u_south + u_north - 4*u_center)*(1/dx)*(1/dx);

	// 			// implement upwinding
	// 			float velocity = get_velocity(x_i);
	// 			if(velocity > 0){
	// 				u_adv = velocity * (u_grid.data[idx] - u_west) / dx;
	// 			} else {
	// 				u_adv = velocity * (u_east - u_grid.data[idx]) / dx;
	// 			}
	// 			// update the concentration data using forward euler
	// 			u_update.data[idx] = u_grid.data[idx] + dt*(diffusion*u_lap + u_adv);
	// 		}
	// 	}
		
	// 	// // write the previous state to a files for later visualization
	// 	// char fpath[25];
	// 	// char* fptr = fpath;
	// 	// sprintf(fptr, "./output/state_%d.csv", n);
	// 	// write_to_file(u_grid, fptr);
		

	// 	int i;
	// 	// overwrite prev iteration with current iteration
	// 	for(i = 0; i <rows*cols; i++){
	// 		u_grid.data[i] = u_update.data[i];
	// 	}

	// }

	free(local_data);
	deallocate_matrix(&u_grid);
	
	MPI_Finalize();

}
