#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "serial.h"

int allocate_vector(Vector* vec, int dim_x, int dim_y){
    // Function to allocate space for a vector
    // Inputs:
    //      vec: vector to allocate space for
    //      dim_x: integer representing number of rows desired
    //      dim_y: integer representing number of columns desired
    // Ouputs:
    //      flag: integer flag, 0 if successful, -1 if not

    // set metadata
    vec->dim_x = dim_x;
    vec->dim_y = dim_y;

    // allocate space
    vec->data = (float*)malloc((vec->dim_x*vec->dim_y)*sizeof(float));

    // ensure allocation was successful
    if(vec->data == NULL){
        fprintf(stderr, "Error allocating data. \n");
        return -1;
    }

    return 0;
}

int create_zero_grid(Vector* vec){
    // Function to fill an initialized vector with zeros
    // Inputs:
    //      vec:  vector to be filled with 0s
    // Outputs:
    //      flag: integer flag, 0 if successful, -1 if not

    // fill the array with 0s
    int i;
    for(i = 0; i < vec->dim_x*vec->dim_y; ++i) vec->data[i] = 0;

    return 0;
}

int create_grid(Vector* vec, float start, float end, char dir){
    // Function to create X and Y grids
    // Inputs:
    //      vec:   vector to fill with the X or Y grid
    //      start: float representing where the grid should start
    //      end:   float representing where the grid should end
    //      dir:   character representing whether it should be an X or Y grid
    // Outputs:
    //      flag:  integer flag, 0 if successful, -1 if not

    // if it should be an X grid, do one thing
    if(dir == 'X' || dir == 'x'){
        // determine the spacing of the grid
        float step = (end - start) / vec->dim_y;

        // loop through the rows
        int row, start_idx, end_idx;
        for(row = 0; row < vec->dim_x; ++row){
            // determine start and end indices of this row
            start_idx = row*vec->dim_y;
            end_idx = (row + 1)*vec->dim_y - 1;

            // declare a float to keep track of current values
            float current_value = start + step / 2;

            // fill in the row
            int idx;
            for(idx = start_idx; idx <= end_idx; ++idx){
                // assign value to current index and update it
                vec->data[idx] = current_value;
                current_value += step;
            }
        }

    }

    // if it should be a Y grid, do another
    else if(dir == 'Y' || dir == 'y'){
        // determine the spacing of the grid
        float step = (end - start) / vec->dim_x;

        // decare a float to keep track of current values
        float current_value = start + step / 2;

        // loop through the rows
        int row, start_idx, end_idx;
        for(row = 0; row < vec->dim_x; ++row){
            // determine start and end indices of this row
            start_idx = row*vec->dim_y;
            end_idx = (row + 1)*vec->dim_y - 1;

            // set the first and last entries of the row
            vec->data[start_idx] = start + step / 2;
            vec->data[end_idx] = end - step / 2;

            // fill in the row
            int idx;
            for(idx = start_idx; idx <= end_idx; ++idx){
                // assign current value to each index
                vec->data[idx] = current_value;
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

Vector square_matrix(Vector* vec){
    // Function to square a vector elementwise
    // Inputs:
    //      vec: vector to square elementwise
    // Outputs:
    //      res: squared vector

    // allocate result vector
    Vector res;
    int flag = allocate_vector(&res, vec->dim_x, vec->dim_y);

    // calculate entries of result vector
    int i;
    for(i = 0; i < res.dim_x*res.dim_y; ++i){
        res.data[i] = vec->data[i]*vec->data[i];
    }

    return res;
}

int subtract_matrices(Vector* vec1, Vector vec2){
    // Function to subtract two matrices and store the result in the first
    // Inputs:
    //      vec1: vector to subtract the second array from
    //      vec2: vector to be subtracted from the first array
    // Outputs:
    //      flag: integer flag, 0 if successful, -1 if not

    // ensure dimensions allign
    if(vec1->dim_x != vec2.dim_x || vec1->dim_y != vec2.dim_y){
        fprintf(stderr, "Error: matrix dimensions must be the same. \n");
        return -1;
    }

    // do the subtraction
    int i;
    for(i = 0; i < vec1->dim_x*vec1->dim_y; ++i) vec1->data[i] -= vec2.data[i];

    return 0;

}

int scalar_multiply(Vector* vec, float c){
    // Function to multiply each element of an array by a scalar
    // Inputs:
    //      vec:  vector whose elements we want to multiply by a scalar
    //      c:    the scalar to multiply by
    // Outputs:
    //      flag: integer flag, 0 if successful, -1 if not

    // do the multiplication
    int i;
    for(i = 0; i < vec->dim_x*vec->dim_y; ++i) vec->data[i] *= c;

    return 0;
}

int initialize_concentration_vector(Vector* u_grid, Vector* x_grid, Vector* y_grid){
    // Function to set a concentration vector to its initial condition
    // Inputs:
    //      u_grid: vector to hold concentration data
    //      x_grid: vector containing an x grid
    //      y_grid: vector contianing a y grid
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

float find_max(Vector vec){
    // Function to find the maximum value in a vector
    // Inputs:
    //      vec: vector to find the maximum value of
    // Outputs:
    //      max: the maximum value encountered in the vector

    float max = -9999.99;
    int i;
    for(i = 0; i < vec.dim_x*vec.dim_y; ++i){
        if(vec.data[i] > max){
            max = vec.data[i];
        }
    }

    return max;
}

float get_velocity(int x){
    // Function to get the velocity corresponsding to a certain x coordinate
    // Currently just a constant function
    // Inputs:
    //      x:   an integer representing where in the grid we are
    // Outputs:
    //      vel: a float representing velocity

    return 0.1;
}

int write_to_file(Vector*){
    // Function to write the data of a matrix to a file
    // Not implemented yet
    
    return 0
}

int deallocate_vector(Vector* vec){

    free(vec->data);
    vec->data = NULL;

    vec->dim_x = 0;
    vec->dim_y = 0;

    return 0;
}
