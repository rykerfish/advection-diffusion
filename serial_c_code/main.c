#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "serial.h"

int main() {

	// declare vectors
	Vector u_grid; 
	Vector x_grid;
	Vector y_grid;
	
	// allocate vectors
	int rows = 128;
	int cols = 128;
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
	
	// set time step parameters
	float cfl = 0.01;
	float dt = cfl * dx;
	float t_final = 10;
	int n_steps = t_final / dt;

	// set physical parameters
	float diffusion = 0.01;

	// do the simulation
	Vector u_update;
	allocate_vector(&u_update, rows, cols);

	// loop through the timesteps
	int n;
	for(n = 0; n < n_steps; ++n){
		float u_lap, u_adv;

		int x_i, y_i;
		for(x_i = 0; x_i < cols; ++x_i){
			for(y_i = 0; y_i < rows; ++y_i){
				int idx = x_i + cols*y_i;
				float u_center = u_grid.data[idx];
				
				// 2D finite laplacian with finite boundary conditions
				float u_west, u_east, u_north, u_south;
				
				if(x_i == 0){ // left edge
					u_west = u_grid.data[idx + cols - 1];
				} else {
					u_west = u_grid.data[idx - 1];
				}

				if(x_i == cols - 1){ // right edge
					u_east = u_grid.data[idx - cols + 1];
				} else {
					u_east = u_grid.data[idx + 1];
				}

				if(y_i == 0){ // top edge
					u_north = u_grid.data[(rows - 1)*cols + x_i];
				} else {
					u_north = u_grid.data[idx - cols];
				}

				if(y_i == rows - 1){ // bottom edge
					u_south = u_grid.data[x_i];
				} else {
					u_south = u_grid.data[idx + cols];
				}

				// calculate finite difference laplacian with a 5 point stencil
				u_lap = (u_west + u_east + u_south + u_north - 4*u_center)*(1/dx)*(1/dx);

				// implement upwinding
				float velocity = get_velocity(x_i);
				if(velocity > 0){
					u_adv = velocity * (u_grid.data[idx] - u_west) / dx;
				} else {
					u_adv = velocity * (u_east - u_grid.data[idx]) / dx;
				}
				// update the concentration data using forward euler
				u_update.data[idx] = u_grid.data[idx] + dt*(diffusion*u_lap + u_adv);
			}
		}
		int i;
		// overwrite prev iteration with current iteration
		for(i = 0; i <rows*cols; i++){
			u_grid.data[i] = u_update.data[i];
		}

	}

	char* filepath = "../out/serial_sim_out.txt";
	write_to_file(u_grid, filepath);
	
	// clean up
	deallocate_vector(&u_grid);
	deallocate_vector(&x_grid);
	deallocate_vector(&y_grid);
	deallocate_vector(&u_update);

	return 0;
}
