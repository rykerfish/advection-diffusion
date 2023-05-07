#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parallel.h"

int main(int argc, char** argv) {

	// declare an integer to be used as a loop variable throughout the program
	int i;

	// ensure we have enough command line arguments
	if(argc < 1){
		printf("Error: Should call with <grid size>\n");
		return 1;
	}

	// interpret command line argument
	int grid_size = atoi(argv[1]);

	// setup our MPI environment and determine number of processors used
	MPI_Init(&argc, &argv);
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// determine current rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Setup initial concentration vector on rank 0, but first declare variables that need to be accessible 
	// by other ranks or outside the if statement
	Matrix u_grid; 
	float dx;
	if(rank == 0){
		// declare matrices
		Matrix x_grid;
		Matrix y_grid;

		// allocate matrices
		int rows = grid_size;
		int cols = grid_size;
		allocate_matrix(&u_grid, cols, rows);
		allocate_matrix(&x_grid, cols, rows);
		allocate_matrix(&y_grid, cols, rows);

		// create grids for constructing the concentration grid
		float a = -3.0; 
		float b = 3.0;
		float x_bounds[] = {a, b};
		float y_bounds[] = {a, b};
		create_grid(&x_grid, x_bounds[0], x_bounds[1], 'x');
		create_grid(&y_grid, y_bounds[0], y_bounds[1], 'y');

		// caluclate grid spacing
		dx = x_grid.data[1] - x_grid.data[0];

		// construct the concentration grid
		initialize_concentration_matrix(&u_grid, &x_grid, &y_grid);

		// clean up
		deallocate_matrix(&x_grid);	
		deallocate_matrix(&y_grid);
	}

	// broadcast dx to all other ranks
	MPI_Bcast(&dx, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// calculate size of subgrids
	int row_len = grid_size;
	int num_rows_in_block = grid_size / nprocs;
	int local_size = row_len * num_rows_in_block; 
	int padded_size = local_size + 2*row_len; // add 2 for ghost region rows

	// allocated memory for local padded array
	float* padded_data = (float*)malloc(padded_size * sizeof(float));
	
	// scatter the data to each process
	// start row_len into local_data so that there is room for ghost region
	MPI_Scatter(u_grid.data, local_size, MPI_FLOAT, &padded_data[row_len], local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// determine neighbors of each rank for use in ghost region communication
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
	float t_final = 50;
	int n_steps = t_final / dt;

	// set physical parameters
	float diffusion = 0.01;
	float velocity = 0.1;

	// declare variable for use in the loop
	int n;
	int plot_step = 20000;
	int img_counter = 0;
	int cols = row_len;
	int rows = num_rows_in_block;

	// declare a matrix for use in updating concentration values
	Matrix u_update;
	allocate_matrix(&u_update, rows, cols);

	// loop through each time step
	for(n = 0; n < n_steps; ++n){
		// ensure all threads make it to the loop
		MPI_Barrier(MPI_COMM_WORLD);

		// send necessary ghost region data
		perform_ghost_comms(padded_data, local_size, row_len, up_neighbor, down_neighbor);

		// declare useful parameters for the inner loops
		int x_i, y_i;
		float u_lap, u_adv;

		// loop over every location in this thread
		for(x_i = 0; x_i < cols; ++x_i){
			for(y_i = 0; y_i < rows; ++y_i){

				// calculate the index of the current location in u_update (idx) and in padded_data (pad_idx)
				int idx = x_i + cols*y_i;
				int pad_idx = idx + cols;
	
				// declare u_west and u_east so we can get them from u_lap
				float u_west, u_east;
				
				// compute laplacian and advection components
				// add 1 to y_i to get past the padded row on top
				u_lap = compute_laplacian(padded_data, x_i, y_i + 1, rows, cols, dx, &u_east, &u_west);
				u_adv = compute_advection(padded_data, pad_idx, velocity, dx, u_east, u_west);

				// store the updated concentration data in u_update using forward euler
				u_update.data[idx] = padded_data[pad_idx] + dt*(diffusion*u_lap + u_adv);
			}
		}

		// overwrite prev iteration with current iteration
		for(i = 0; i <rows*cols; i++){
			padded_data[i+cols] = u_update.data[i];
		}

		//  write the previous state to a files for later visualization
		if(n % plot_step == 0){
			MPI_Gather(&padded_data[row_len], local_size, MPI_FLOAT, u_grid.data, local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
			if(rank == 0){
				char fpath[25];
				char* fptr = fpath;
				sprintf(fptr, "./output/state_%d.csv", img_counter);
				img_counter++;
				write_to_file(u_grid, fptr);
			}
		}
	}

	// combine all of the data onto thread 0 to be written out
	MPI_Gather(&padded_data[row_len], local_size, MPI_FLOAT, u_grid.data, local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// write the final concentration matrix to a txt file
	if(rank == 0){
		char* filepath = "../test_out/parallel_sim_out.txt";
		write_to_file(u_grid, filepath);

		deallocate_matrix(&u_grid);
	}

	// clean up
	free(padded_data);
	deallocate_matrix(&u_update);
	MPI_Finalize();

}
