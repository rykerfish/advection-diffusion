#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "serial.h"

int main() {

	Vector u_grid; 
	Vector x_grid;
	Vector y_grid;
	
	allocate_vector(&u_grid, 100, 100);
	allocate_vector(&x_grid, 100, 100);
	allocate_vector(&y_grid, 100, 100);
	
	create_grid(&x_grid, -3, 3, 'x');
	create_grid(&y_grid, -3, 3, 'y');

	initialize_concentration_vector(&u_grid, &x_grid, &y_grid);
	
	int i;
	for(i = 0; i < 100; ++i){
		int j;
		for(j = 0; j < 100; ++j){
			printf("%f ", (&u_grid)->data[j*100 + i]);
		}
		printf("\n");
	}

	return 0;

}
