#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parallel.h"

int main(int argc, char** argv) {

	int i;

	if(argc < 1){
		printf("Error: Should call with <grid size>\n");
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
	float dx;
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

		dx = x_grid.data[1] - x_grid.data[0];

		// construct the concentration grid
		initialize_concentration_matrix(&u_grid, &x_grid, &y_grid);

		deallocate_matrix(&x_grid);	
		deallocate_matrix(&y_grid);
	}

	MPI_Bcast(&dx, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// calculate size of subgrid
	int row_len = grid_size;
	int num_rows_in_block = grid_size / nprocs;
	int local_size = row_len * num_rows_in_block; 
	int padded_size = local_size + 2*row_len; // add 2 for ghost region rows

	// allocated memory for local array
	float* padded_data = (float*)malloc(padded_size * sizeof(float));
	
	// scatter the data to each process
	// start row_len into local_data so that there is room for ghost region
	MPI_Scatter(u_grid.data, local_size, MPI_FLOAT, &padded_data[row_len], local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	int up_neighbor, down_neighbor;
	if(rank == 0){
		up_neighbor = nprocs - 1;
	} else {
		up_neighbor = rank - 1;
	}

	if(rank == nprocs - 1){
		down_neighbor = 0;
	} else {
		down_neighbor = rank + 1;
	}

	// set time step parameters
	float cfl = 0.01;
	float dt = cfl * dx;
	float t_final = 10;
	int n_steps = t_final / dt;

	// set physical parameters
	float diffusion = 0.01;
	float velocity = 0.1;

	// do the simulation
	Matrix u_update;

	// loop through the timesteps
	int n;
	int plot_step = 1000;
	int cols = row_len;
	int rows = num_rows_in_block;

	printf("n steps: %d\n", n_steps);

	allocate_matrix(&u_update, rows, cols);

	for(n = 0; n < n_steps; ++n){
		MPI_Barrier(MPI_COMM_WORLD);

		float u_lap, u_adv;
		perform_ghost_comms(padded_data, local_size, row_len, up_neighbor, down_neighbor);

		int x_i, y_i;

		for(x_i = 0; x_i < cols; ++x_i){
			for(y_i = 0; y_i < rows; ++y_i){
				int idx = x_i + cols*y_i;
				int pad_idx = idx + cols;
	
				float u_west, u_east;
				
				u_lap = compute_laplacian(padded_data, x_i, y_i + 1, rows, cols, dx, &u_east, &u_west);

				if(rank == 0){
					if(x_i == 0 && y_i == 0){
						printf("0,0 u_lap, %f\n", u_lap);
					} else if(x_i == 25 && y_i == 21){
						printf("25,21 u_lap, %f\n", u_lap);
					}
				}


				// implement upwinding
				if(velocity > 0){
					u_adv = velocity * (padded_data[pad_idx] - u_west) / dx;
				} else {
					u_adv = velocity * (u_east - padded_data[pad_idx]) / dx;
				}
				// update the concentration data using forward euler
				u_update.data[idx] = padded_data[pad_idx] + dt*(diffusion*u_lap + u_adv);
			}
		}

		// overwrite prev iteration with current iteration
		for(i = 0; i <rows*cols; i++){
			padded_data[i+cols] = u_update.data[i];
		}

		// write the previous state to a files for later visualization
		// if(n % plot_step == 0){
		// 	char fpath[25];
		// 	char* fptr = fpath;
		// 	sprintf(fptr, "./output/state_%d.csv", n);
		// 	write_to_file(u_grid, fptr);
		// }
	}
	

	MPI_Gather(&padded_data[row_len], local_size, MPI_FLOAT, u_grid.data, local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	if(rank == 0){
		char* filepath = "../out/parallel_sim_out.txt";
		write_to_file(u_grid, filepath);

		deallocate_matrix(&u_grid);
	}


	free(padded_data);
	deallocate_matrix(&u_update);

	
	MPI_Finalize();

}
