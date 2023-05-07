# Loose files
- test.sh: the main test file that should get submitted to slurm. Calls the other test files in the parallel/ and serial/ directories.
- integration_tests.py: has the integration tests. These were written in python for ease of use as files are read in and compared, which is more painful in C.
- *_scaling.sh: The requested scaling test files.

# Parallel
- main.c: main simulation file, written in parallel.
- parallel.*: .c and .h files containing the function implementations used in the simulation.
- test_parallel.c: test file containing tests written in parallel.
- test.sh: bash script to submit to slurm to test parallel code.
- output: Directory holding data files to compute visualizations.

# Serial
- main.c: main simulation file, written in serial.
- serial.*: .c and .h files containing the function implementations used in the simulation.
- test.c: test file containing tests written in serial.
- test.sh: bash script to submit to slurm to test serial code.

# Slurm_out
- No files at present. Gets filled with output files when code/test.sh is submitted to slurm. 

# Test_out
- No files at present. Gets filled with output files when code/test.sh is submitted to slurm. 

# Matlab
- advection.m: original advection-diffusion simulation, used to check correctness of serial code. Not used in project, but included for completeness.
- test_sim.m and test.m: tests used to check correctness of serial code implementation. Not used in project, but included for completeness.