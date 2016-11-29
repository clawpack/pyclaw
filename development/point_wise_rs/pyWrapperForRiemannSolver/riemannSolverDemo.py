
from __future__ import absolute_import
import RiemannSolver
from RiemannSolver import *
from pylab import plot,  figure, suptitle

rs = RiemannSolver(timeSteps =1, mwaves = 2, mx = 800, meqn = 2, maux = 2)
    
rs.q = random.random((rs.meqn,rs.mx))
rs.aux = random.random((rs.maux,rs.mx))
waves1, s1  = rs.solveVectorized(timer = True)
waves2, s2  = rs.solvePointwize(timer = True)



figure(1)
waveIndex = 0
componentIndex = 0

suptitle("Figure 1: shows the curve of wave number {0}, for the component number {1}".format(waveIndex+1, componentIndex+2))	
plot(waves1[componentIndex,waveIndex,:], "g")
