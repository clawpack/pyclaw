import numpy as np
import pylab as pl

#Script to solve 1D advection equation:
#
# q_t + q_x = 0
#
#Using first-order finite differences.

#Grid:
m=100
x=np.linspace(0,1,m+1); x=x[:-1]
h=x[1]-x[0]
cflnum=0.95
k=cflnum*h
T=2.
N=round(T/k)

# Initial condition:
u = np.exp(-10*(x-0.5)**2)

for n in xrange(N+1):
    u[1:]=u[1:]-cflnum*(u[1:]-u[:-1])
    #Uncomment below if you want it to plot:
    #pl.clf(); pl.hold(False)
    #pl.plot(x,u); pl.axis([-0,1,-0.1,1.1]); pl.draw()
