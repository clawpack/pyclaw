import RiemannSolver
from RiemannSolver import *
from numpy import load, absolute


tolerance = 0.00005
print("This script asserts the results for input read from data files q.npy and aux.npy with tolerance =", tolerance)



rs = RiemannSolver( 10, 2, 2**15, 2, 2)
    
rs.q = load("q.npy")
rs.aux = load("aux.npy")

wavesRead = load("waves.npy")
sRead = load("s.npy")

waves1, s1  = rs.solveVectorized(timer = True)
#waves2, s2  = rs.solvePointwize(timer = True)

assert (absolute(wavesRead - waves1)< tolerance).all()
assert (absolute(sRead -s1)< tolerance).all()
