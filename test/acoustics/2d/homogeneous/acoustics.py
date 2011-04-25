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


def acoustics2D(use_PETSc=True,kernel_language='Fortran',iplot=False,petscPlot=False,useController=True,htmlplot=False):
    """
    Example python script for solving the 2d acoustics equations.
    """
    if use_PETSc:
        from petsc4py import PETSc
        import petclaw as myclaw
        output_format='petsc'
        from petclaw.evolve.clawpack import PetClawSolver2D as mySolver
    else: #Pure pyclaw
        import pyclaw as myclaw
        output_format='ascii'
        from pyclaw.evolve.clawpack import ClawSolver2D as mySolver

    from pyclaw.solution import Solution
    from pyclaw.controller import Controller

    from petclaw import plot

    # Initialize grid
    mx=100; my=100
    x = myclaw.grid.Dimension('x',-1.0,1.0,mx,mthbc_lower=1,mthbc_upper=1)
    y = myclaw.grid.Dimension('y',-1.0,1.0,my,mthbc_lower=1,mthbc_upper=1)
    grid = myclaw.grid.Grid([x,y])

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
    tfinal = 0.12
    qinit(grid)
    initial_solution = Solution(grid)

    solver = mySolver()
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 2
    solver.mthlim = [4]*solver.mwaves
    solver.dt=np.min(grid.d)/grid.aux_global['cc']*solver.cfl_desired

    claw = Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.output_format = 'petsc'
    claw.tfinal = tfinal
    claw.solutions['n'] = initial_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    if htmlplot:  plot.plotHTML()
    if petscPlot: plot.plotPetsc(claw)
    if iplot:     plot.plotInteractive()

    pressure=claw.frames[claw.nout].grid.gqVec.getArray().reshape([grid.local_n[0],grid.local_n[1],grid.meqn])[:,:,0]
    return pressure


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        acoustics2D(*args,**kwargs)
    else: acoustics2D()
