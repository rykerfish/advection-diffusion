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
	int rows = 100;
	int cols = 100;
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

	// construct the concentration grid
	initialize_concentration_vector(&u_grid, &x_grid, &y_grid);
	
	// print grid to console for comparison
	int i;
	for(i = 0; i < 100; ++i){
		int j;
		for(j = 0; j < 100; ++j){
			printf("%f ", (&u_grid)->data[j*100 + i]);
		}
		printf("\n");
	}
	
	// clean up
	deallocate_vector(&u_grid);
	deallocate_vector(&x_grid);
	deallocate_vector(&y_grid);
	
	return 0;

}
