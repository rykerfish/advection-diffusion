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

    printf("PADDED: %d\n", padded_size);

	// allocated memory for local array
	local_data = (float*)malloc(padded_size * sizeof(float));
	
	// scatter the data to each process
	// start row_len into local_data so that there is room for ghost region
	MPI_Scatter(u_grid, local_size, MPI_FLOAT, local_data + row_len, local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    printf("data: %f", local_data[3]);

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

int deallocate_matrix(Matrix* mat){

    free(mat->data);
    mat->data = NULL;

    mat->dim_x = 0;
    mat->dim_y = 0;

    return 0;
}
