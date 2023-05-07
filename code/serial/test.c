#include <stdio.h>
#include <stdlib.h>
#include "serial.h"

int main(){

    unsigned int testCount = 0;
	unsigned int testsPassed = 0;
	unsigned int testsFailed = 0;


    // ------------ Testing matrix creation/teardown --------
    printf("--------- BEGIN SERIAL UNIT TESTS -----------\n");
    testCount++;
    printf("Test %d. Testing matrix creation/teardown\n", testCount);
    
    int grid_size = 5;
    Matrix temp_mat;
    int flag = allocate_matrix(&temp_mat, grid_size, grid_size);
    if(flag != 0){
        printf("WARNING: Test %d failed due to bad matrix allocation.\n", testCount);
        testsFailed++;
    } else if(temp_mat.dim_x != grid_size || temp_mat.dim_y != grid_size){
        printf("WARNING: Test %d failed due to incorrect matrix dimensions.\n", testCount);
        testsFailed++;
    } else {
        // check that entries can be set. WARNING: Could cause failure
        temp_mat.data[0] = 10;
        temp_mat.data[grid_size*grid_size-1] = 0.1;

        deallocate_matrix(&temp_mat);

        if(temp_mat.data != NULL || temp_mat.dim_x != 0 || temp_mat.dim_y != 0){
            printf("WARNING: Test %d failed due to bad matrix deallocation.\n", testCount);
            testsFailed++;
        } else{
            testsPassed++;
            printf("Test %d passed.\n", testCount);
        }
    }

    // ------------ Test X and Y grid setup -----------

    testCount++;
    printf("Test %d. Testing X and Y grid creation\n", testCount);

    grid_size = 200;
    unsigned int test_failed = 0;
    Matrix x_grid, y_grid;
    allocate_matrix(&x_grid, grid_size, grid_size);
    allocate_matrix(&y_grid, grid_size, grid_size);
    create_grid(&x_grid, -3, 3, 'x');
    create_grid(&y_grid, -3, 3, 'y');

    int i;
    for(i = 0; i < grid_size*grid_size; ++i){
        if(test_failed != 0) break;
        // test x grid values
        if(i % grid_size == 0){
            if(abs(x_grid.data[i] + 2.985) > 1e-6) test_failed = 1;
        } else if(i % grid_size == 62){
            if(abs(x_grid.data[i] + 1.125) > 1e-6) test_failed = 1;
        } else if(i % grid_size == 183){
            if(abs(x_grid.data[i] - 2.505) > 1e-6) test_failed = 1;
        }
        // test y grid values
        if(i < grid_size){
            if(abs(y_grid.data[i] + 2.985) > 1e-6) test_failed = 2;
        } else if(i >= grid_size*62 && i < grid_size*63){
            if(abs(y_grid.data[i] + 1.125) > 1e-6) test_failed = 2;
        } else if(i >= grid_size*183 && i < grid_size*184){
            if(abs(y_grid.data[i] - 2.505) > 1e-6) test_failed = 2;
        }
    }

    if(test_failed != 0){
        testsFailed++;
        if(test_failed == 1){
            printf("WARNING: Test %d failed due to incorrect x grid setup.\n", testCount);
        } else if(test_failed == 2){
            printf("WARNING: Test %d failed due to incorrect y grid setup.\n", testCount);
        }
    } else{
        testsPassed++;
        printf("Test %d passed.\n", testCount);
    }

    printf("------------ SERIAL UNIT TESTS RESULTS-------------\n");
    printf("Tests passed: %d out of %d\n", testsPassed, testCount);
    printf("Tests failed: %d out of %d\n", testsFailed, testCount);

    // Use output of X and Y grids to test concentration grid construction
    test_failed = 0;
    Matrix u_grid;
    allocate_matrix(&u_grid, grid_size, grid_size);
    initialize_concentration_matrix(&u_grid, &x_grid, &y_grid);

    if(abs(u_grid.data[104+186*grid_size] - 0.003504744309636) > 1e-6) test_failed = 1;
    printf("--------- BEGIN SERIAL INTEGRATION TESTS -----------\n");

    printf("Creating serial output files for integration tests:\n");

    int rows = 6;
    int cols = 5;
    Matrix integration_mat;
    allocate_matrix(&integration_mat, cols, rows);
    int x_i, y_i;
    for(y_i = 0; y_i < rows; y_i++){
        for(x_i = 0; x_i < cols; x_i++){
            integration_mat.data[x_i + y_i*cols] = x_i + y_i*cols;
        }
    }

    Matrix u_lap;
    Matrix u_adv;
    allocate_matrix(&u_lap, cols, rows);
    allocate_matrix(&u_adv, cols, rows);
    float east_temp;
    float west_temp;
    for(y_i = 0; y_i < rows; y_i++){
        for(x_i = 0; x_i < cols; x_i++){
            u_lap.data[x_i + y_i*cols] = compute_laplacian(integration_mat.data, x_i, y_i, rows, cols, 1, &east_temp, &west_temp);
            u_adv.data[x_i + y_i*cols] = compute_advection(integration_mat.data, x_i + y_i*cols, 1, 1, east_temp, west_temp);
        }
    }


    write_to_file(u_lap, "../test_out/integration_laplacian.txt");
    printf("Serial laplacian for integration tests written to test_out/integration_laplacian.txt\n");
    write_to_file(u_adv, "../test_out/integration_advection.txt");
    printf("Serial advection for integration tests written to test_out/integration_advection.txt\n");

    deallocate_matrix(&integration_mat);
    deallocate_matrix(&u_adv);
    deallocate_matrix(&u_lap);


    return 0;
}