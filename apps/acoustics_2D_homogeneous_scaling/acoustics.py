#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,width=0.2):
    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    q[1,:,:] = 0.
    q[2,:,:] = 0.
    grid.q=q


def acoustics2D(iplot=False,petscPlot=False,useController=True,htmlplot=False):
    """
    Example python script for solving the 2d acoustics equations.
    """

    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.clawpack import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

    # Initialize grid
    # Amal: need to modify mx, my final time to make fixed number of time steps
    size = PETSc.Comm.getSize(PETSc.COMM_WORLD)
    mx=4*int(np.sqrt(size*10000)); my=mx
    rank = PETSc.Comm.getRank(PETSc.COMM_WORLD)
    if rank == 0:
        print "mx, my = ",mx, my
    x = Dimension('x',-1.0,1.0,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',-1.0,1.0,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])

    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    grid.aux_global['rho']= rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']= zz
    grid.aux_global['cc']=cc
    from dimsp2 import cparam
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    grid.meqn = 3
    grid.mbc = 2
    tfinal = 0.2/np.sqrt(size)
    qinit(grid)
    initial_solution = Solution(grid)

    solver = PetClawSolver2D()
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 2
    solver.mthlim = [4]*solver.mwaves
    solver.dt=np.min(grid.d)/grid.aux_global['cc']*solver.cfl_desired
    sol = {"n":initial_solution}
    solver.setup(sol)
    solver.evolve_to_time(sol,tfinal/10.) 
     
    claw = Controller()
    claw.nout = 1
    claw.keep_copy = True
    size = PETSc.Comm.getSize(PETSc.COMM_WORLD)
    claw.outdir = './_output_'+str(size)+'/'
    # The output format MUST be set to petsc!
    claw.output_format = 'petsc'
    claw.tfinal = tfinal
    claw.solutions['n'] = initial_solution
    claw.solver = solver

    # Solve
    import time
    start=time.time()    
    status = claw.run()
    end=time.time()
    duration1 = end-start
    rank = PETSc.Comm.getRank(PETSc.COMM_WORLD)
    print 'controller.run took'+str(duration1)+' seconds, for process '+str(rank)
    if rank ==0:
        print 'number of steps: '+ str(claw.solver.status.get('numsteps'))
                            
if __name__=="__main__":
    import sys
    from petsc4py import PETSc
    WithArgs = False
    generateProfile = False
    proccessesList = [0,5]
    
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        WithArgs= True

    if generateProfile:
       rank =PETSc.Comm.getRank(PETSc.COMM_WORLD)
       size =PETSc.Comm.getSize(PETSc.COMM_WORLD)
       if rank in proccessesList:
           import cProfile
           if WithArgs: cProfile.run('acoustics2D(*args,**kwargs)', 'profile'+str(rank)+'_'+str(size))
           else: cProfile.run('acoustics2D()', 'profile'+str(rank)+'_'+str(size))
       else:
           print "process"+str(rank) +"not profiled"
           if WithArgs: acoustics2D(*args,**kwargs)
           else: acoustics2D()
    else:
        if WithArgs: acoustics2D(*args,**kwargs)
        else: acoustics2D()
                   
