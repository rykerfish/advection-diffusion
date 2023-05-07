import numpy as np

print("-------RUNNING INTEGRATION TESTS----------")
testCount = 0
testsPassed = 0
testsFailed = 0

testCount = testCount + 1

par = np.genfromtxt('./test_out/integration_laplacian_parallel.txt', dtype=float, delimiter=',')
serial = np.genfromtxt('./test_out/integration_laplacian.txt', dtype=float, delimiter=',')

norm = sum(sum(abs(par - serial)))
if(norm < 1e-4):
    testsPassed = testsPassed + 1
    print("Laplacian comparison passed: parallel results sufficiently close to serial code.")
else:
    print("WARNING: Full simulation test failed. Norm too big: ", norm)
    testsFailed = testsFailed + 1

testCount = testCount + 1

par = np.genfromtxt('./test_out/integration_advection_parallel.txt', dtype=float, delimiter=',')
serial = np.genfromtxt('./test_out/integration_advection.txt', dtype=float, delimiter=',')

norm = sum(sum(abs(par - serial)))
if(norm < 1e-4):
    testsPassed = testsPassed + 1
    print("Laplacian comparison passed: parallel results sufficiently close to serial code.")
else:
    print("WARNING: Full simulation test failed. Norm too big: ", norm)
    testsFailed = testsFailed + 1

testCount = testCount + 1

par = np.genfromtxt('./test_out/parallel_sim_out.txt', dtype=float, delimiter=',')
serial = np.genfromtxt('./test_out/serial_sim_out.txt', dtype=float, delimiter=',')

norm = sum(sum(abs(par - serial)))
if(norm < 1e-4):
    testsPassed = testsPassed + 1
    print("Full simulation test passed: parallel results sufficiently close to serial code.")
else:
    print("WARNING: Full simulation test failed. Norm too big: ", norm)
    testsFailed = testsFailed + 1


print("TEST RESULTS")
print("Tests passed: " + str(testsPassed) + " out of " + str(testCount))
print("Tests failed: " + str(testsFailed) + " out of " + str(testCount))