#ifndef __SERIAL_H__
#define __SERIAL_H__

typedef struct Vector_tag {
    float* data; // vector to hold concentration (or other) data
    int dim_x; // number of rows
    int dim_y; // number of columns
} Vector;

// Function to allocate space for a vector
// Inputs:
//      vec:   vector to allocate space for
//      dim_x: integer representing number of rows desired
//      dim_y: integer representing number of columns desired
// Ouputs:
//      flag:  integer flag, 0 if successful, -1 if not
int allocate_vector(Vector* vec, int dim_x, int dim_y);

// Function to fill an initialized vector with zeros
// Inputs:
//      vec:  vector to be filled with 0s
// Outputs:
//      flag: integer flag, 0 if successful, -1 if not
int create_zero_grid(Vector* vec);

// Function to create X and Y grids
// Inputs:
//      vec:   vector to fill with the X or Y grid
//      start: float representing where the grid should start
//      end:   float representing where the grid should end
//      dir:   character representing whether it should be an X or Y grid
// Outputs:
//      flag:  integer flag, 0 if successful, -1 if not
int create_grid(Vector* vec, float start, float end, char dir);

// Function to square a vector elementwise
// Inputs:
//      vec: vector to square elementwise
// Outputs:
//      res: squared vector
Vector square_matrix(Vector* vec);

// Function to subtract two matrices and store the result in the first
// Inputs:
//      vec1: vector to subtract the second array from
//      vec2: vector to be subtracted from the first array
// Outputs:
//      flag: integer flag, 0 if successful, -1 if not
int subtract_matrices(Vector* vec1, Vector vec2);

// Function to multiply each element of an array by a scalar
// Inputs:
//      vec:  vector whose elements we want to multiply by a scalar
//      c:    the scalar to multiply by
// Outputs:
//      flag: integer flag, 0 if successful, -1 if not
int scalar_multiply(Vector* vec, float c);

// Function to set a concentration vector to its initial condition
// Inputs:
//      u_grid: vector to hold concentration data
//      x_grid: vector containing an x grid
//      y_grid: vector contianing a y grid
// Outputs:
//      flag:  integer flag, 0 if successful, -1 if not
int initialize_concentration_vector(Vector* u_grid, Vector* x_grid, Vector* y_grid);

// Function to find the maximum value in a vector
// Inputs:
//      vec: vector to find the maximum value of
// Outputs:
//      max: the maximum value encountered in the vector
float find_max(Vector vec);

// Function to get the velocity corresponsding to a certain x coordinate
// Currently just a constant function
// Inputs:
//      x:   an integer representing where in the grid we are
// Outputs:
//      vel: a float representing velocity
float get_velocity(int x);

// Function to write the data of a matrix to a file
// Inputs:
//      vec:      vector whose contents should be written
//      filepath: string representing path to file
// Outputs:
//      flag:     integer flag, 0 if successful, -1 otherwise
int write_to_file(Vector vec, char* filepath);

// Function to cleanup memory from a vector
// Inputs:
//      vec: vector to be deallocated.
// Outputs:
//      None, but vec will be set to null
int deallocate_vector(Vector* vec);

#endif
