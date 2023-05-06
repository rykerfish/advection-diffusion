#include <stdio.h>
#include <stdlib.h>
#include "parallel.h"

int main(){

    unsigned int testCount = 0;
	unsigned int testsPassed = 0;
	unsigned int testsFailed = 0;


    // ------------ Testing matrix creation/teardown --------
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

    // ------------ Testing grid setup -----------








    return 0;
}