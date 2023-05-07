import numpy as np

par = np.genfromtxt('./test_out/integration_laplacian_parallel.txt', dtype=float, delimiter=',')
serial = np.genfromtxt('./test_out/integration_laplacian.txt', dtype=float, delimiter=',')

norm = sum(sum(abs(par - serial)))
if(norm < 1e-4):
    print("Laplacian comparison passed: parallel results sufficiently close to serial code.")
else:
    print("WARNING: Full simulation test failed. Norm too big: ", norm)

par = np.genfromtxt('./test_out/integration_advection_parallel.txt', dtype=float, delimiter=',')
serial = np.genfromtxt('./test_out/integration_advection.txt', dtype=float, delimiter=',')

norm = sum(sum(abs(par - serial)))
if(norm < 1e-4):
    print("Laplacian comparison passed: parallel results sufficiently close to serial code.")
else:
    print("WARNING: Full simulation test failed. Norm too big: ", norm)

par = np.genfromtxt('./test_out/parallel_sim_out.txt', dtype=float, delimiter=',')
serial = np.genfromtxt('./test_out/serial_sim_out.txt', dtype=float, delimiter=',')

norm = sum(sum(abs(par - serial)))
if(norm < 1e-4):
    print("Full simulation test passed: parallel results sufficiently close to serial code.")
else:
    print("WARNING: Full simulation test failed. Norm too big: ", norm)