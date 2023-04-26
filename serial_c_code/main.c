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
	int n_rows = 5;
	int n_cols = 5;
	allocate_vector(&u_grid, n_rows, n_cols);
	allocate_vector(&x_grid, n_rows, n_cols);
	allocate_vector(&y_grid, n_rows, n_cols);
	
	// create grids for constructing the concentration grid
	float x_start = -3.0;
	float x_end = 3.0;
	float y_start = -3.0;
	float y_end = 3.0;
	create_grid(&x_grid, x_start, x_end, 'x');
	create_grid(&y_grid, y_start, y_end, 'y');

	float dx = x_grid.data[1] - x_grid.data[0];
	float dy = y_grid.data[n_rows + 1] - y_grid.data[0];

	// construct the concentration grid
	initialize_concentration_vector(&u_grid, &x_grid, &y_grid);
	
	// print grid to console for comparison
	int i;
	// for(i = 0; i < 100; ++i){
	// 	int j;
	// 	for(j = 0; j < 100; ++j){
	// 		printf("%f ", (&u_grid)->data[j*100 + i]);
	// 	}
	// 	printf("\n");
	// }
	// return 0;

	// time stepping parameters
	float cfl = 0.01;
	float dt = cfl*dx;
	float T_final = 10;
	int N_steps = T_final/dt;

	// physical parameters
	float vel = 0.1; // velocity
	float D = 0.01; // diffusion

	int n;


	float u_center, u_west, u_east, u_north, u_south;

	Vector u_update;
	allocate_vector(&u_update, n_rows, n_cols);

	for(n=0; n<N_steps; n++){
		int y, x;
		float u_lap, u_adv;
		for(y=0; y < n_rows; y++){
			for(x=0; x < n_cols; x++){

				// LAPLACIAN PART

				u_center = u_grid.data[x + y*n_cols];

				if(x == 0){ // left edge
					u_west = u_grid.data[(n_cols-1) + y*n_cols];
				} else{
					u_west = u_grid.data[y*n_cols + x - 1];
				} 

				if(x == n_cols-1){ // right edge
					u_east = u_grid.data[0 + y*n_cols];
				} else{
					u_east = u_grid.data[(x + 1) + y*n_cols];
				}

				if(y == 0){ // top edge
					u_north = u_grid.data[x + (n_rows-1)*n_cols];
				} else{
					u_north = u_grid.data[x + (y-1)*n_cols];
				}

				if(y == n_rows-1){ // bottom edge
					u_south = u_grid.data[x + 0];
				} else{
					u_south = u_grid.data[x + (y+1)*n_cols];
				}

				u_lap = (u_west + u_east + u_south + u_north - 4*u_center)*(1/(dx*dx));

				// UPWINDING PART
				if(vel > 0){
					u_adv = vel*(u_center - u_west)/dx;
				} else{
					u_adv = vel*(u_east - u_center)/dx;
				}
				printf("%f ", u_adv);
			}
			printf("\n");

			// forward euler for single cell
			u_update.data[x + y*n_cols] = u_grid.data[x + y*n_cols] + dt*D*u_lap + dt*u_adv;
		}
		return 0;

		for(i = 0; i < n_rows*n_cols; i++){
			u_grid.data[i] = u_update.data[i];
		}

		// for(y = 0; y < n_rows; ++y){
		// 	for(x = 0; x < n_cols; ++x){
		// 		printf("%f ", u_update.data[x*n_cols + y]);
		// 	}
		// 	printf("\n");
		// }
	

	}
	int x, y;

	
	// clean up
	deallocate_vector(&u_grid);
	deallocate_vector(&x_grid);
	deallocate_vector(&y_grid);
	deallocate_vector(&u_update);
	
	return 0;

}
