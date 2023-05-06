#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <mpi.h>
#include "parallel.h"

int allocate_matrix(Matrix* mat, int dim_x, int dim_y){
    // Function to allocate space for a matrix
    // Inputs:
    //      mat: matrix to allocate space for
    //      dim_x: integer representing number of rows desired
    //      dim_y: integer representing number of columns desired
    // Ouputs:
    //      flag: integer flag, 0 if successful, -1 if not

    // set metadata
    mat->dim_x = dim_x;
    mat->dim_y = dim_y;

    // allocate space
    mat->data = (float*)malloc((mat->dim_x*mat->dim_y)*sizeof(float));

    // ensure allocation was successful
    if(mat->data == NULL){
        fprintf(stderr, "Error allocating data. \n");
        return -1;
    }

    return 0;
}

int create_grid(Matrix* mat, float start, float end, char dir){
    // Function to create X and Y grids
    // Inputs:
    //      mat:   matrix to fill with the X or Y grid
    //      start: float representing where the grid should start
    //      end:   float representing where the grid should end
    //      dir:   character representing whether it should be an X or Y grid
    // Outputs:
    //      flag:  integer flag, 0 if successful, -1 if not

    // if it should be an X grid, do one thing
    if(dir == 'X' || dir == 'x'){
        // determine the spacing of the grid
        float step = (end - start) / mat->dim_y;

        // loop through the rows
        int row, start_idx, end_idx;
        for(row = 0; row < mat->dim_x; ++row){
            // determine start and end indices of this row
            start_idx = row*mat->dim_y;
            end_idx = (row + 1)*mat->dim_y - 1;

            // declare a float to keep track of current values
            float current_value = start + step / 2;

            // fill in the row
            int idx;
            for(idx = start_idx; idx <= end_idx; ++idx){
                // assign value to current index and update it
                mat->data[idx] = current_value;
                current_value += step;
            }
        }

    }

    // if it should be a Y grid, do another
    else if(dir == 'Y' || dir == 'y'){
        // determine the spacing of the grid
        float step = (end - start) / mat->dim_x;

        // decare a float to keep track of current values
        float current_value = start + step / 2;

        // loop through the rows
        int row, start_idx, end_idx;
        for(row = 0; row < mat->dim_x; ++row){
            // determine start and end indices of this row
            start_idx = row*mat->dim_y;
            end_idx = (row + 1)*mat->dim_y - 1;

            // set the first and last entries of the row
            mat->data[start_idx] = start + step / 2;
            mat->data[end_idx] = end - step / 2;

            // fill in the row
            int idx;
            for(idx = start_idx; idx <= end_idx; ++idx){
                // assign current value to each index
                mat->data[idx] = current_value;
            }

            // update current value
            current_value += step;
        }
    }

    // otherwise, throw an error
    else{
        fprintf(stderr, "Error: must enter either X/x or Y/y. \n");
        return 0;
    }

    return 1;
}

Matrix square_matrix(Matrix* mat){
    // Function to square a matrix elementwise
    // Inputs:
    //      mat: matrix to square elementwise
    // Outputs:
    //      res: squared matrix

    // allocate result matrix
    Matrix res;
    int flag = allocate_matrix(&res, mat->dim_x, mat->dim_y);

    // calculate entries of result matrix
    int i;
    for(i = 0; i < res.dim_x*res.dim_y; ++i){
        res.data[i] = mat->data[i]*mat->data[i];
    }

    return res;
}

int subtract_matrices(Matrix* mat1, Matrix mat2){
    // Function to subtract two matrices and store the result in the first
    // Inputs:
    //      mat1: matrix to subtract the second array from
    //      mat2: matrix to be subtracted from the first array
    // Outputs:
    //      flag: integer flag, 0 if successful, -1 if not

    // ensure dimensions allign
    if(mat1->dim_x != mat2.dim_x || mat1->dim_y != mat2.dim_y){
        fprintf(stderr, "Error: matrix dimensions must be the same. \n");
        return -1;
    }

    // do the subtraction
    int i;
    for(i = 0; i < mat1->dim_x*mat1->dim_y; ++i) mat1->data[i] -= mat2.data[i];

    return 0;

}


int initialize_concentration_matrix(Matrix* u_grid, Matrix* x_grid, Matrix* y_grid){
    // Function to set a concentration matrix to its initial condition
    // Inputs:
    //      u_grid: matrix to hold concentration data
    //      x_grid: matrix containing an x grid
    //      y_grid: matrix contianing a y grid
    // Outputs:
    //      flag:  integer flag, 0 if successful, -1 if not

    // make sure dimensions allign
    if(u_grid->dim_x != x_grid->dim_x || x_grid->dim_x != y_grid->dim_x || u_grid->dim_y != x_grid->dim_y || x_grid->dim_y != y_grid->dim_x){
        fprintf(stderr, "Error: dimensions do not allign. \n");
        return -1;
    }

    // set all data of concentration data to 0
    int i;
    for(i = 0; i < u_grid->dim_x*u_grid->dim_y; ++i) u_grid->data[i] = 0;

    // use x and y grids to construct the initial condition
    int flag;
    flag = subtract_matrices(u_grid, square_matrix(x_grid));
    flag = subtract_matrices(u_grid, square_matrix(y_grid));

    for(i = 0; i < u_grid->dim_x*u_grid->dim_y; ++i) u_grid->data[i] = 3*expf(u_grid->data[i]);

    return 0;

}

int perform_scatter(float* u_grid, int grid_size, int nprocs, float* local_data){

    // calculate size of subgrid
	int row_len = grid_size;
	int num_rows_in_block = grid_size / nprocs;
	int local_size = row_len * num_rows_in_block; 
	int padded_size = local_size + 2*row_len; // add 2 for ghost region rows

	// allocate memory for local array
	local_data = (float*)malloc(padded_size * sizeof(float));
	
	// scatter the data to each process
	// start row_len into local_data so that there is room for ghost region
	MPI_Scatter(u_grid, local_size, MPI_FLOAT, &local_data[row_len], local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return local_size;
}

void perform_ghost_comms(float* local_data, int local_size, int row_len, int up_neighbor, int down_neighbor){
    
    // setup for communcations
	MPI_Status status[4];
	MPI_Request request[4];
	int i;
	for(i = 0; i < 4; i++) request[i] = MPI_REQUEST_NULL;

    int padded_size = local_size + 2*row_len;

	// do ghost region comms
	MPI_Isend(local_data + row_len, row_len, MPI_FLOAT, up_neighbor, 0, MPI_COMM_WORLD, &request[0]); // send first row of data to upper neighbor
	MPI_Isend(local_data + local_size, row_len, MPI_FLOAT, down_neighbor, 0, MPI_COMM_WORLD, &request[1]); // send last row of data to lower neighbor
	MPI_Irecv(local_data, row_len, MPI_FLOAT, up_neighbor, 0, MPI_COMM_WORLD, &request[2]); // recv data from upper neighbor
	MPI_Irecv(local_data + padded_size - row_len, row_len, MPI_FLOAT, down_neighbor, 0, MPI_COMM_WORLD, &request[3]); // recv data from lower neighbor

	MPI_Waitall(4, request, status);
}

int write_to_file(Matrix mat, char* filepath){
    // Function to write the data of a matrix to a file
    // Inputs:
    //      mat:      matrix whose contents should be written
    //      filepath: string representing path to file
    // Outputs:
    //      flag:     integer flag, 0 if successful, -1 otherwise

    // declare file pointer
    FILE* fptr;

    // open the file and ensure it was successful
    fptr = fopen(filepath, "w");
    if(fptr == NULL){
        fprintf(stderr, "Error opening file, error number %d \n", errno);
        return -1;
    }

    int i, j;
    for(i = 0; i < mat.dim_x; ++i){
        for(j = 0; j < mat.dim_y; ++j){
            fprintf(fptr, "%f", mat.data[i + j*mat.dim_x]);
            if(j < mat.dim_y - 1) fprintf(fptr, ", ");
        }
        fprintf(fptr, "\n");
    }
    
    fclose(fptr);

    return 0;
}

float compute_laplacian(float* data, int x_i, int y_i, int rows, int cols, float dx, float* u_east, float* u_west){

    int idx = x_i + cols * y_i;
    float u_north, u_south, u_center;
    int i;
    u_center = data[idx];


    if(x_i == 0){ // left edge
		*u_west = data[idx + cols - 1];
	} else {
		*u_west = data[idx - 1];
	}

	if(x_i == cols - 1){ // right edge
		*u_east = data[idx - cols + 1];
	} else {
		*u_east = data[idx + 1];
	}

	if(y_i == 0){ // top edge
		u_north = data[(rows - 1)*cols + x_i];
	} else {
		u_north = data[idx - cols];
	}

	if(y_i == rows - 1){ // bottom edge
		u_south = data[x_i];
	} else {
		u_south = data[idx + cols];
	}

	// calculate finite difference laplacian with a 5 point stencil
    float u_lap = (*u_west + *u_east + u_south + u_north - 4.0*u_center)*(1.0/dx)*(1.0/dx);
    return u_lap;
}

float compute_advection(float* data, int idx, float velocity, float dx, float u_east, float u_west){
    
    float u_adv;

    if(velocity > 0){
        u_adv = velocity * (data[idx] - u_west) / dx;
    } else {
        u_adv = velocity * (u_east - data[idx]) / dx;
    }

    return u_adv;
}

int deallocate_matrix(Matrix* mat){

    free(mat->data);
    mat->data = NULL;

    mat->dim_x = 0;
    mat->dim_y = 0;

    return 0;
}
