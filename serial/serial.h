#ifndef __SERIAL_H__
#define __SERIAL_H__

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

// Function to scatter data onto each thread
//
// Inputs:
//      u_grid:     pointer to an array of floats to be scattered (assumed to be on rank 0)
//      grid_size:  integer representing width and height of the array (must be square)
//      nprocs:     integer representing the number of processes
//      local_data: pointer to an array on each thread to store the scattered data in
// Ouptus:
//      local_size: integer representing the number of floats given to each thread
int perform_scatter(float* u_grid, int grid_size, int nprocs, float* local_data);

// Function to perform communications between threads to fill in the ghost regions
// 
// Inputs:
//      local_data:    pointer to a padded array on each thread to send/recieve information
//      local_size:    integer represnting the number of floats in local_data, not counting the padding
//      row_len:       integer representing the length of each row in local_data
//      up_neighbor:   integer representing which thread is "above" this thread
//      down_neighbor: integer representing which thread is "below" this thread
// Outputs:
//      None
void perform_ghost_comms(float* local_data, int local_size, int row_len, int up_neighbor, int down_neighbor);

// Function to write the data of a matrix to a file
// Inputs:
//      vec:      matrix whose contents should be written
//      filepath: string representing path to file
// Outputs:
//      flag:     integer flag, 0 if successful, -1 otherwise
int write_to_file(Matrix mat, char* filepath);

// Function to compute the location of a point
// Inputs:
//      data:   pointer to an array of floats containing the point to compute the laplacian of
//      x_i:    the x coordinate (column) of the point within data whose laplacian should be calculated
//      y_i:    the y coordinate (row) of the point within data whose laplacian should be calculated
//      rows:   the number of rows in data (used for boundary conditions)
//      cols:   the number of columns in data (used for boundary conditions)
//      dx:     the spacing between points in data
//      u_east: a pointer to a float in which the value of the point to the east will be stored
//      u_west: a pointer to a float in which the value of the point to the west will be stored
// Outputs:
//      u_lap:  float representing the laplacian at the current point
float compute_laplacian(float* data, int x_i, int y_i, int rows, int cols, float dx, float* u_east, float* u_west);

// Function to compute the advection at the current point
// Inputs:
//      data:     pointer to an array of floats containing the point to compute the advection at
//      idx:      the index of the current point
//      velocity: the velocity at the current point
//      dx:       the spacing between grid points
//      u_east:   the value of the point to the east of the current point
//      u_west:   the value of the point to the west of the current point
// Outputs:
//      u_adv:    the advection at the current point
float compute_advection(float* data, int idx, float velocity, float dx, float u_east, float u_west);

// Function to deallocate a matrix
// Inputs: 
//      mat:  the matrix to be deallocated
// Outpus:
//      flag: integer flag, 0 if successful, -1 otherwise
int deallocate_matrix(Matrix* mat);

#endif