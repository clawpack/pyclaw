#!/usr/bin/env python
import numpy as np
from petsc4py import PETSc
from six.moves import range

class PetCLAW:
    def advection1D(self, M, cfl, T):
        '''Script to solve 1D advection equation:
        q_t + q_x = 0
        Using first-order finite differences'''
        
        # PETSc DA object, which handles structured grids, here are
        # requesting M global points 
        da = PETSc.DA().create([M])
        print(da.getSizes())

        ranges = da.getRanges()
        print(ranges[0][1] - ranges[0][0])
        

        # this solves the problem on the domain [0 1] by default
        da.setUniformCoordinates() 

        # view commands dump to the screen
        da.view()

        # this gets the coordinate vector, this is akin to the
        # linspace(0,1,M) command 
        xvec = da.getCoordinates()
        xvec.view()
        # access the local data array within the globally distributed
        # vector 
        x = xvec.getArray()

        h = x[1]-x[0]
        k=cfl*h

        # global vector represents coordinated data with other
        # processors 
        fg = da.createGlobalVector()
        # local vector represents local data that is not shared
        fl = da.createLocalVector()

        # we will operate on the local vector, then coordinate with
        # other processors on the global vector 

        N=int(round(T/k))
        
        # Initial condition:
        q = np.exp(-10*(x-0.5)**2)
        # this is a little tricky.  the local vector contains ghost
        # points, which represents data that comes in from the global
        # vector when we are coordinating.  On the other hand, the
        # global vector only ever contains the data we own.  So when
        # we set initial conditions, we do it on the global vector,
        # then scatter to the local  
        fg.setArray(q)
        da.globalToLocal(fg,fl)

        # this should dump the properly set initial conditions
        fg.view()
        
        for n in range(N+1):
            # grab the working array out of the local vector
            # (including ghost points)
            q = fl.getArray()

            # operate on the local data
            q[1:]=q[1:]-cfl*(q[1:]-q[:-1])

            # restore the working array
            fl.setArray(q)

            # this is a local update, the local array in the local
            # vector is copied into the corresponding local array in
            # the global vector, ghost points are discarded 
            da.localToGlobal(fl,fg)
            
            fg.view()
            # this is a global update coordinated over all processes.
            # The local array in the global vector is copied into the
            # corresponding local array in the local vector.
            # In addition, the ghost values needed to operate locally
            # are sent over MPI to the correct positions in the local
            # array in the local vector.
            da.globalToLocal(fg,fl)

    def run(self):
        OptDB = PETSc.Options()
        M = OptDB.getInt('M', 16)
        cfl = OptDB.getReal('cfl',0.95)
        T = OptDB.getReal('T',2.)
        self.advection1D(M, cfl, T)
        print('Done')
        return
          
if __name__ == '__main__':
    PetCLAW().run()
          

    
    
