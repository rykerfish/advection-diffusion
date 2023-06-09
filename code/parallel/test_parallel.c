#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parallel.h"

int main(){


    MPI_Init(NULL, NULL);

    unsigned int testCount = 0;
    unsigned int testsPassed = 0;
    unsigned int testsFailed = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) printf("\n\n---------------BEGIN PARALLEL UNIT TESTS-------------\n");

    // ------------ testing scatter -------------
    Matrix mat;
    int grid_size = 3;
    int i;
    if(rank == 0){
        testCount++;
        printf("Test %d. Testing scatter.\n", testCount);
        allocate_matrix(&mat, grid_size, grid_size);
        
        for(i = 0; i < 9; i++){
            mat.data[i] = i;
        }
    }

    float* local_data = NULL;
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // calculate size of subgrid
    int row_len = grid_size;
    int num_rows_in_block = grid_size / nprocs;
    int local_size = row_len * num_rows_in_block; 
    int padded_size = local_size + 2*row_len; // add 2 for ghost region rows

    // allocate memory for local array
    local_data = (float*)malloc(padded_size * sizeof(float));
	
    // scatter the data to each process
    // start row_len into local_data so that there is room for ghost region
    MPI_Scatter(mat.data, local_size, MPI_FLOAT, &local_data[row_len], local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    int elements_wrong[3] = {0, 0, 0};
    for(i = 0; i < grid_size; i++){
     	float target = 3*rank + i;
    	int local_idx = i + grid_size;
        if(local_data[local_idx] != target){
            elements_wrong[rank] += 1;
        }
    }

    int sum_wrong = 0;
    for(i = 0; i < 3; i++){
        if(elements_wrong[rank] != 0){
            printf("TEST FAILED: Rank %d has the wrong data after scatter\n", rank);
            sum_wrong++;
        }
    }

    if(sum_wrong > 0){
        testsFailed++;
    } else{
        if(rank == 0){
            testsPassed++;
            printf("Test %d passed.\n", testCount);
        }

    }

    if(testsFailed > 1) testsFailed = 1;



    // ------------- Testing ghost comms ---------
    if(rank == 0){
        testCount++;
        printf("Test %d. Testing ghost region communications.\n", testCount);
    }
 
    int up_neighbor, down_neighbor;

    if(rank == 0){
	    up_neighbor = nprocs - 1;
    } else{
	    up_neighbor = rank - 1;
    }

    if(rank == nprocs - 1){
	    down_neighbor = 0;
    } else{
	    down_neighbor = rank + 1;
    }

    perform_ghost_comms(local_data, local_size, grid_size, up_neighbor, down_neighbor);

    int entries_wrong = 0;
    if(rank == 1){
        for(i = 0; i < 9; ++i){
            if(local_data[i] != i) entries_wrong++;
        }
    }

    MPI_Bcast(&entries_wrong, 1, MPI_INT, 1, MPI_COMM_WORLD);

    if(rank == 0){
        if(entries_wrong > 0){
            testsFailed++;
            printf("TEST FAILED: Rank 1 has the wrong data after ghost region communication. \n");
        } 
        else{
            testsPassed++;
            printf("Test %d passed.\n", testCount);
        }
    }

    // ----------- Testing Laplacian in the middle-------------------
    if(rank == 0){
        testCount++;
        printf("Test %d. Testing laplacian calculation in middle of array.\n", testCount);
    }

    int test_passed = 0;
    if(rank == 1){
        float target_lap = 0.0;
        float u_west, u_east;
        float actual_lap = compute_laplacian(local_data, 1, 1, 3, 3, (float)1, &u_east, &u_west);
        if(actual_lap == target_lap) test_passed = 1;
    }

    MPI_Bcast(&test_passed, 1, MPI_INT, 1, MPI_COMM_WORLD);

    if(rank == 0){
        if(test_passed == 0){
            testsFailed++;
            printf("TEST FAILED: Incorrect laplacian using middle of rank 1 \n");
        } else{
            testsPassed++;
            printf("Test %d passed.\n", testCount);
        }
    }

    // ----------- Testing Laplacian on the side-------------------
    if(rank == 0){
        testCount++;
        printf("Test %d. Testing laplacian calculation on the side of array.\n", testCount);
    }

    test_passed = 0;
    if(rank == 1){
        float target_lap = -3.0;
        float u_west, u_east;
        float actual_lap = compute_laplacian(local_data, 2, 1, 3, 3, (float)1, &u_east, &u_west);
        if(actual_lap == target_lap) test_passed = 1;
    }

    MPI_Bcast(&test_passed, 1, MPI_INT, 1, MPI_COMM_WORLD);

    if(rank == 0){
        if(test_passed == 0){
            testsFailed++;
            printf("TEST FAILED: Incorrect laplacian using right edge of rank 1 \n");
        } else{
            testsPassed++;
            printf("Test %d passed.\n", testCount);
        }
    }

    // ----------- Testing Laplacian on the corner-------------------
    if(rank == 0){
        testCount++;
        printf("Test %d. Testing laplacian calculation on the corner of array.\n", testCount);
    }

    test_passed = 0;
    if(rank == 1){
        float target_lap = 6.0;
        float u_west, u_east;
        float actual_lap = compute_laplacian(local_data, 2, 0, 3, 3, (float)1, &u_east, &u_west);
        if(actual_lap == target_lap) test_passed = 1;
    }

    MPI_Bcast(&test_passed, 1, MPI_INT, 1, MPI_COMM_WORLD);

    if(rank == 0){
        if(test_passed == 0){
            testsFailed++;
            printf("TEST FAILED: Incorrect laplacian using top right corner of rank 1 \n");
        } else{
            testsPassed++;
            printf("Test %d passed.\n", testCount);
        }
    }

    // ------------ Testing advective term in middle ------------
    
    MPI_Barrier(MPI_COMM_WORLD);
    testCount++;
    if(rank == 0) printf("Test %d. Testing advection term calculation.\n", testCount);

    float target_adv = 1.0;
    int x_i = 1; 
    int y_i = 0;
    int idx = x_i + grid_size*(y_i + 1);
    float u_west = 3;
    float u_east = 5;

    
    MPI_Barrier(MPI_COMM_WORLD);


    test_passed = 0;
    if(rank == 1){
        float actual_adv = compute_advection(local_data, idx, 1, 1, u_east, u_west);

        if(abs(actual_adv-target_adv) < 1e-5){
            test_passed++;
        }
    }
    MPI_Bcast(&test_passed, 1, MPI_INT, 1, MPI_COMM_WORLD);

    if(rank == 0){
        if(test_passed != 1){
            printf("WARNING: Test %d failed: incorrect advective term.\n", testCount);
            testsFailed++;
        } else {
            testsPassed++;
            printf("Test %d passed.\n", testCount);
        }
    }

    if(rank == 0){
        printf("------------ PARALLEL UNIT TESTS RESULTS-------------\n");
        printf("Tests passed: %d out of %d\n", testsPassed, testCount);
        printf("Tests failed: %d out of %d\n", testsFailed, testCount);

        printf("--------- BEGIN PARALLEL INTEGRATION TESTS -----------\n");
        printf("Creating parallel output files for integration tests:\n");
    } 

    int rows = 6;
    int cols = 5;
    Matrix integration_mat;
    if(rank == 0){
        allocate_matrix(&integration_mat, cols, rows);
        for(y_i = 0; y_i < rows; y_i++){
            for(x_i = 0; x_i < cols; x_i++){
                integration_mat.data[x_i + y_i*cols] = x_i + y_i*cols;
            }
        }
    }

    // calculate size of subgrid
    num_rows_in_block = rows / nprocs;
    local_size = cols * num_rows_in_block;
    padded_size = local_size + 2*cols; // add 2 for ghost region rows
	
    // allocate memory for local array
    float* integration_local = (float*)malloc(padded_size * sizeof(float));

	if(rank == 0){
		up_neighbor = nprocs - 1;
	} else {
		up_neighbor = rank - 1;
	}

	if(rank == nprocs - 1){
		down_neighbor = 0;
	} else {
		down_neighbor = rank + 1;
	}

    // scatter the data to each process
    // start row_len into local_data so that there is room for ghost region
    MPI_Scatter(integration_mat.data, local_size, MPI_FLOAT, &integration_local[cols], local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    perform_ghost_comms(integration_local, local_size, cols, up_neighbor, down_neighbor);

    Matrix u_lap;
    Matrix u_adv;
    if(rank == 0){
        allocate_matrix(&u_lap, cols, rows);
        allocate_matrix(&u_adv, cols, rows);
    }

    local_size = cols*num_rows_in_block;

    float* lap_local = malloc(local_size * sizeof(float));
    float* adv_local = malloc(local_size * sizeof(float));

    float east_temp, west_temp;
    for(y_i = 0; y_i < num_rows_in_block; y_i++){
        for(x_i = 0; x_i < cols; x_i++){
            lap_local[x_i + y_i*cols] = compute_laplacian(integration_local, x_i, y_i+1, num_rows_in_block, cols, 1, &east_temp, &west_temp);
            adv_local[x_i + y_i*cols] = compute_advection(integration_local, x_i+(y_i+1)*cols, 1, 1, east_temp, west_temp);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(lap_local, local_size, MPI_FLOAT, u_lap.data, local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(adv_local, local_size, MPI_FLOAT, u_adv.data, local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    free(integration_local);
    free(lap_local);
    free(adv_local);

    if(rank == 0){
        write_to_file(u_lap, "../test_out/integration_laplacian_parallel.txt");
        printf("Parallel laplacian for integration tests written to in test_out/integration_laplacian_parallel.txt\n");

        write_to_file(u_adv, "../test_out/integration_advection_parallel.txt");
        printf("Parallel advection for integration tests written to test_out/integration_advection_parallel.txt\n");
        deallocate_matrix(&u_lap);
        deallocate_matrix(&integration_mat);
    }







    if(rank == 0){

	    printf("Done with testing, %d out of %d passed. \n", testsPassed, testCount);

    }

    

    // free(local_data);

    MPI_Finalize();

}
