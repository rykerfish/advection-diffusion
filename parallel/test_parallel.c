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
        printf("------------------------------------------------- \n");
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
        testsPassed++;
    }

    if(testsFailed > 1) testsFailed = 1;



    // ------------- Testing ghost comms ---------
    if(rank == 0){
        testCount++;
        printf("Test %d. Testing ghost region communications.\n", testCount);
        printf("------------------------------------------------- \n");
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
        }
    }


    if(rank == 0){

	    printf("Done with testing, %d out of %d passed. \n", testsPassed, testCount);

    }

    MPI_Finalize();

}
