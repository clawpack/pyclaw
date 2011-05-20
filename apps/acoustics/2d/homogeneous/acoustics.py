#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,width=0.2):
    
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


def acoustics2D(iplot=False,htmlplot=False,use_PETSc=False,outdir='./_output',soltype='classic'):
    """
    Example python script for solving the 2d acoustics equations.
    """

    if use_PETSc:
        from petclaw.grid import Dimension, Grid
        if soltype=='classic':
            from petclaw.evolve.clawpack import PetClawSolver2D as mySolver
        elif soltype=='sharpclaw':
            from petclaw.evolve.sharpclaw import PetSharpClawSolver2D as mySolver
        from petclaw.controller import Controller
    else:
        from pyclaw.grid import Dimension, Grid
        if soltype=='classic':
            from pyclaw.evolve.clawpack import ClawSolver2D as mySolver
        elif soltype=='sharpclaw':
            from pyclaw.evolve.sharpclaw import SharpClawSolver2D as mySolver
        from pyclaw.controller import Controller

    from pyclaw.solution import Solution
    from petclaw import plot

    solver = mySolver()
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 2
    solver.mthlim = [4]*solver.mwaves


    # Initialize grid
    mx=100; my=100
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
    if soltype=='classic':
        from dimsp2 import cparam
    elif soltype=='sharpclaw':
        from flux2 import cparam
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    solver.dt=np.min(grid.d)/grid.aux_global['cc']*solver.cfl_desired

    grid.meqn = 3
    grid.mbc = solver.mbc
    tfinal = 0.12

    if use_PETSc:
        # Initialize petsc Structures for q
        grid.init_q_petsc_structures()

    qinit(grid)
    initial_solution = Solution(grid)

    claw = Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.tfinal = tfinal
    claw.solutions['n'] = initial_solution
    claw.solver = solver
    claw.outdir=outdir

    # Solve
    status = claw.run()

    if htmlplot:  plot.plotHTML(outdir=outdir,format=claw.output_format)
    if iplot:     plot.plotInteractive(outdir=outdir,format=claw.output_format)

    if use_PETSc:
        pressure=claw.frames[claw.nout].grid.gqVec.getArray().reshape([grid.local_n[0],grid.local_n[1],grid.meqn])[:,:,0]
    else:
        pressure=claw.frames[claw.nout].grid.q[:,:,0]
    return pressure


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        acoustics2D(*args,**kwargs)
    else: acoustics2D()
