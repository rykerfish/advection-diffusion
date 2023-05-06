#ifndef __PARALLEL_H__
#define __PARALLEL_H__

typedef struct Matrix_tag {
    float* data; // matrix to hold concentration (or other) data
    int dim_x; // number of rows
    int dim_y; // number of columns
} Matrix;

// Function to allocate space for a matrix
// Inputs:
//      mat:   matrix to allocate space for
//      dim_x: integer representing number of rows desired
//      dim_y: integer representing number of columns desired
// Ouputs:
//      flag:  integer flag, 0 if successful, -1 if not
int allocate_matrix(Matrix* mat, int dim_x, int dim_y);

// Function to create X and Y grids
// Inputs:
//      mat:   matrix to fill with the X or Y grid
//      start: float representing where the grid should start
//      end:   float representing where the grid should end
//      dir:   character representing whether it should be an X or Y grid
// Outputs:
//      flag:  integer flag, 0 if successful, -1 if not
int create_grid(Matrix* mat, float start, float end, char dir);

// Function to square a matrix elementwise
// Inputs:
//      mat: matrix to square elementwise
// Outputs:
//      res: squared matrix
Matrix square_matrix(Matrix* mat);

// Function to subtract two matrices and store the result in the first
// Inputs:
//      mat1: matrix to subtract the second array from
//      mat2: matrix to be subtracted from the first array
// Outputs:
//      flag: integer flag, 0 if successful, -1 if not
int subtract_matrices(Matrix* mat1, Matrix mat2);


// Function to set a concentration matrix to its initial condition
// Inputs:
//      u_grid: matrix to hold concentration data
//      x_grid: matrix containing an x grid
//      y_grid: matrix contianing a y grid
// Outputs:
//      flag:  integer flag, 0 if successful, -1 if not
int initialize_concentration_matrix(Matrix* u_grid, Matrix* x_grid, Matrix* y_grid);

void perform_ghost_comms(float* local_data, int local_size, int row_len, int up_neighbor, int down_neighbor);

int perform_scatter(float* u_grid, int grid_size, int nprocs, float* local_data);

// Function to write the data of a matrix to a file
// Inputs:
//      vec:      matrix whose contents should be written
//      filepath: string representing path to file
// Outputs:
//      flag:     integer flag, 0 if successful, -1 otherwise
int write_to_file(Matrix mat, char* filepath);

float compute_laplacian(float* data, int x_i, int y_i, int rows, int cols, float dx, float* u_east, float* u_west);

float calculate_advection(float* data, int idx, float velocity, float dx, float u_west, float u_east);

// Function to cleanup memory from a matrix
// Inputs:
//      mat: matrix to be deallocated.
// Outputs:
//      None, but mat will be set to null
int deallocate_matrix(Matrix* mat);

#endif
