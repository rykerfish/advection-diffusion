## Advection and Diffusion Equation Solution
### Ryker Fish and Daniel Callan Richards-Dinger

### Writeup Questions (answers must be complete sentences)


1. (5%) Describe the serial algorithm that you are parallelizing. What problem is it aiming to address? Describe the algorithm with specific sentences or pseudocode. How many operations does the serial algorithm for a particular problem size?

    a. We chose to parallelize a solution to the advection diffusion equation that uses upwinding for the advective component and the laplacian for the diffusive component. The advection diffusion equations describes a physical situation in which particles (or any other physical quantity) are moved due to two processes: diffusion and advection. For example, it can be used to model the spread of methane gas due to wind (advection) and natural diffusion. Psuedocode for our solution is shown below:

    + Set up an evenely spaced grid of initial concentrations (we chose a 2d gaussian)
    + Loop through the time steps. At each time step:
        + Estimate the laplacian at each grid location using a 5 point stencil with periodic boundary conditions
        + 

2. (10%) Using detailed sentences or pseudocode, how are you breaking up this algorithm using a distributed memory (MPI) model? Why is this algorithm not embarrassingly parallel? Describe the initial distribution of the data/model, and what data need to be transferred between processes throughout. If it is helpful to introduce new notation, please do so, but be sure to define each index/variable specifically defined at its first use. You may include diagrams as needed. 

3. (5%) Describe your test suite and use its results to make a case that your code is behaving as expected. How do you know your code is right? Are there any limitations to your test suite? 

4. (6%) Describe the aspect of your code you chose to visualize to understand (and quickly check) whether the code is be behaving as expected. What did you visualize? How did you do this? How do you know it should look like this? Include some of the figures you produced. 

5. (7%) Describe how you designed your strong and weak scalability tests (dimensions and number of processes). Fill out two tables with the headers (number of procs, weak test time, weak test efficiency) and (number of procs, strong test time, strong test efficiency). Comment on whether your efficiency changes, and make a case for what the biggest contributing factors may be to any slowdowns. 
