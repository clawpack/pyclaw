from numpy import random, empty,zeros
#from numpy import *
from FRS import vec_rp, pw_rp
from time import clock
from pylab import plot,  figure
from six.moves import range



class RiemannSolver:
  

    def __init__(self, timeSteps, mwaves, mx, meqn, maux, q = None, aux = None):
       
        self.timeSteps = timeSteps
        self.mwaves = mwaves
        self.mx = mx
        self.meqn = meqn
        self.maux = maux
        if q is not None:
            self.q = q
        if aux is not None:
            self.aux = aux

    def solveVectorized(self, timer = None):
       
        if timer is not None:
            print("Solving... vectorized")
            start = clock()

        for counter in range(self.timeSteps):
            waves, s = vec_rp(self.mwaves,self.q[:,0:self.mx-1],self.q[:,1:self.mx],self.aux[:,0:self.mx-1],self.aux[:,1:self.mx])

        if timer is not None:
            end = clock()
            self.elapsedTime = end-start
            print("elapsed time: ", self.elapsedTime, " seconds")

        return waves, s




    def solvePointwize(self, timer = None):

        timeSteps, mwaves, mx, meqn, maux, q , aux = self.timeSteps, self.mwaves, self.mx, self.meqn, self.maux, self.q , self.aux
        waves = zeros((meqn, mwaves, mx))
        s = zeros((mwaves, mx))
        if timer is not None:
            print("Solving... point-wize")
            start = clock()

        for counter in range(timeSteps):
            for i in range(1,mx-1):  #from 1 to mx-1 unlike fortran (from 1 to mx-2 means including mx-2)
                waves[:,:,i], s[:,i] = pw_rp(mwaves,q[:,i],q[:,i+1],aux[:,i],aux[:,i+1]) # instead of pw_rp(mwaves,q[:,i-1],q[:,i],aux[:,i-1],aux[:,i]) due to the
                                                                                         # differences in indexing among fortran and python

        if timer is not None:
            end = clock()
            self.elapsedTime = end-start
            print("elapsed time: ", self.elapsedTime, " seconds")

        return waves, s
        

        

if __name__ == "__main__":
    rs = RiemannSolver( 1, 2, 2**10, 2, 2)
    
    rs.q = random.random((rs.meqn,rs.mx))   
    rs.aux = random.random((rs.maux,rs.mx))

    
    waves1, s1  = rs.solveVectorized(timer = True)

    
    waves2, s2  = rs.solvePointwize(timer = True)
   
