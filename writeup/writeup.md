## Advection and Diffusion Equation Solution
### Ryker Fish and Daniel Callan Richards-Dinger

### Writeup Questions (answers must be complete sentences)


1. (5%) Describe the serial algorithm that you are parallelizing. What problem is it aiming to address? Describe the algorithm with specific sentences or pseudocode. How many operations does the serial algorithm for a particular problem size?

    a. We chose to parallelize a solution to the advection diffusion equation that uses upwinding for the advective component and the laplacian for the diffusive component. The advection diffusion equations describes a physical situation in which particles (or any other physical quantity) are moved due to two processes: diffusion and advection. For example, it can be used to model the spread of methane gas due to wind (advection) and natural diffusion. Psuedocode for our solution is shown below:

    + Set up an evenely spaced grid of initial concentrations (we chose a 2d gaussian)
    + Loop through the time steps. At each time step:
        + Estimate the laplacian at each grid location using a 5 point stencil with periodic boundary conditions
        + Calculate the advective component using upwinding (multiplying the velocity by the difference between each concentration and the concentration of the point upwind of it)
        + Combine these to get an updated concentration grid (new concentration = old + time step * diffusion * laplcian + time step * advection)

    The number of operations performed for the serial algorithm depends on both grid size and time steps. For each time step, 16 operations are performed at each grid location. Assuming a square grid (which out code does), this results in 16\*t\*n\*n operations, where t is the number of time steps and n is the size of the grid.

2. (10%) Using detailed sentences or pseudocode, how are you breaking up this algorithm using a distributed memory (MPI) model? Why is this algorithm not embarrassingly parallel? Describe the initial distribution of the data/model, and what data need to be transferred between processes throughout. If it is helpful to introduce new notation, please do so, but be sure to define each index/variable specifically defined at its first use. You may include diagrams as needed. 

    a. Our simulation occurs on a grid (say of size N x N), so our paralellization strategy involves decomosing this domain onto multiple threads. Assume that we have np threads. That means that each thread gets an N/np x N matrix to work with. Since the laplacian uses a 5-point stencil, each thread also needs one row of padding above and below it, so each thread's domain was padded to (N/np + 2) x N. Then within the for loop, each thread computes the diffusive and advective component of its data in parallel, and between each time step the ghost regions are updated with new information from the other threads. Finally, whenever the full concentration data is needed (e.g. for visualization), the data from all the threads is gathered and the full matrix can be written to a file.

3. (5%) Describe your test suite and use its results to make a case that your code is behaving as expected. How do you know your code is right? Are there any limitations to your test suite? 

    a. Our test suite is broken up into two main parts: parallel tests and serial tests. The serial tests ensure that the setup of the code works as expected. In this section, we test matrix allocation and deallocation and the setup of our grids, which are used to create the initial concentration vector. We use the output of the grid test to test the setup of the concentration vector, forming our first integration test. All of these grid tests are performed using a spot check method, where enough elements are compared to the known correct values (taken from Matlab) that there will be an extrememly low likelihood of them all being correct if the entire matrix is not. In this section of the test suite we also compute some matrices in serial to use for our integration tests in the parallel section.

    Our parallel section starts with unit tests to test our parallel functions. We begin by testing our scatter function, using a small matrix to ensure the right data ends up on each thread. WE then test our ghost region commnication function, using the same small matrix to ensure the correct padding is used. We then move on to testing our laplacian function, running three unit tests to ensure it computes the laplacian correctly in the middle, on the side, and on a corner of a small test matrix whose laplacians we compute by hand to check. We also test the laplacian on an entire, larger, grid. Our final unit test ensures that our advection works as intended, verified by calculating the expected result by hand. We then perform an integration test on the laplacian, creating a small test matrix, scattering it, adding ghost regions, and then finally computing the laplacian and comparing it to the serial result calculated in the previous test section. 

4. (6%) Describe the aspect of your code you chose to visualize to understand (and quickly check) whether the code is be behaving as expected. What did you visualize? How did you do this? How do you know it should look like this? Include some of the figures you produced. 

    a. Since our code generates concentration grids, we chose to visualize these grids at various time steps as heat maps using Python. We compared these to similar images generated from the serial simulation to verify that they were correct. 

5. (7%) Describe how you designed your strong and weak scalability tests (dimensions and number of processes). Fill out two tables with the headers (number of procs, weak test time, weak test efficiency) and (number of procs, strong test time, strong test efficiency). Comment on whether your efficiency changes, and make a case for what the biggest contributing factors may be to any slowdowns. 

    a. We designed our strong scaling test by first choosing a grid size of 256 x 256 somewhat arbitrarily (although we chose it to be apower of 2). We then chose number of processes from 1 to 16, following the powers of 2. For the weak scaling test, since the operations completed are O(n^2), we chose numbers of processes as 1, 4, and 16 (going up by powers of 4), and grid sizes of 64, 128, and 256, increasing by powers of 2. 

    | # Threads | Strong Test Time | Strong Test Efficiency |
    |-----------|------------------|------------------------|
    | 1         | 2m13.847s        | 1                      |
    | 2         | 56.762s          | 1.179                  |
    | 4         | 27.883s          | 1.2                    |
    | 8         | 14.264s          | 1.1729                 |
    | 16        | 8.063s           | 1.0375                 |

    | # Threads | Weak Test Time | Weak Test Efficiency |
    |-----------|----------------|----------------------|
    | 1         | 1.787s         | 1                    |
    | 4         | 3.178s         | 0.5623               |
    | 16        | 7.262s         | 0.246                |

    Interestingly, in our strong scaling test we see efficiencies of over 1 for all numbers of processes. While this is somewhat surprising, it is not unheard of -- we got similar results in a lab. On the other hand, our efficiency decreased significantly as the number of threads increased in our weak scaling test.
