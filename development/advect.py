#!/usr/bin/env python
import sys
from six.moves import range

try:
    import numpy as np
    from petsc4py import PETSc
except:
    sys.path.append("/opt/share/ksl/petsc4py/dev-aug29/ppc450d/lib/python/")
    sys.path.append("/opt/share/ksl/numpy/dev-aug29/ppc450d/lib/python/")
    
    import numpy as np
    from petsc4py import PETSc

class PetCLAW:
    def advection1D(self, M, cfl, T):
        '''Script to solve 1D advection equation:
        q_t + q_x = 0
        Using first-order finite differences'''
        
        da = PETSc.DA().create([M])
        da.setUniformCoordinates() # solves the problem from 0 to 1
        da.view()

        xvec = da.getCoordinates()
        xvec.view()
        x = xvec.getArray()
        
        h = x[1]-x[0]
        k=cfl*h
        
        fg = da.createGlobalVector()
        fl = da.createLocalVector()
        
        N=int(round(T/k))
        
        # Initial condition:
        q = np.exp(-10*(x-0.5)**2)
        fg.setArray(q)
        da.globalToLocal(fg,fl)

        fg.view()
        
        for n in range(N+1):
            q = fl.getArray()
            q[1:]=q[1:]-cfl*(q[1:]-q[:-1])
            fl.setArray(q)
            da.localToGlobal(fl,fg)            
            fg.view()
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
          

    
    
