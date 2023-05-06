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
    if(rank == 0){
        printf("\n\n---------------PARALLEL TESTS-------------\n\n");
    }
    


    // ------------ testing scatter -------------
    Matrix mat;
    int grid_size = 3;
    int i;
    if(rank == 0){
        testCount++;
        printf("Test %d. Testing scatter\n", testCount);

        allocate_matrix(&mat, grid_size, grid_size);
        
        for(i = 0; i < 9; i++){
            mat.data[i] = i;
            printf("%d\n", mat.data[i]);
        }
        printf("%d\n", mat.data[2]);
        printf("------------------\n");
    }

    float* local_data = NULL;
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // int local_size = perform_scatter(mat.data, grid_size, nprocs, local_data);
    // calculate size of subgrid
    int row_len = grid_size;
    int num_rows_in_block = grid_size / nprocs;
    int local_size = row_len * num_rows_in_block; 
    int padded_size = local_size + 2*row_len; // add 2 for ghost region rows

    // allocate memory for local array
    local_data = (float*)malloc(padded_size * sizeof(float));

    //int local_size;
    //local_size = perform_scatter(mat.data, grid_size, nprocs, local_data);
    //printf("Rank %d done with scatter. \n", rank);	
    // scatter the data to each process
    // start row_len into local_data so that there is room for ghost region
    MPI_Scatter(mat.data, local_size, MPI_FLOAT, &local_data[row_len], local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    int elements_wrong[3] = {0, 0, 0};
    for(i = 0; i < grid_size; i++){
     	float target = 3*rank + i;
    	int local_idx = i + grid_size;
    	fprintf(stderr, "local_idx: %d, data: %f. \n", local_idx, local_data[local_idx]);
        if(local_data[local_idx] != target){
            elements_wrong[rank] += 1;
        }
    }

    if(rank == 2){
        for(i = 0; i < 9; i++){
            fprintf(stderr, "rank 2 value %d, %f\n", i, local_data[i]);
        }
    }

    if(rank == 0){
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
            testsPassed++;
        }
    }


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

    if(rank == 2){
	for(i = 0; i < 9; ++i){
	    fprintf(stderr, "Rank 2 value %d: %f \n", i, local_data[i]); 
	}
    }


    if(rank == 0){

	printf("Done with testing, %d out of %d passed. \n", testsPassed, testCount);

    }

    MPI_Finalize();

}
