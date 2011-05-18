#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,rad=1.0):
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    q[0,:,:] = 0.25*np.pi + 3.25*np.pi*(r<=rad)
    grid.q=q


def kpp(use_PETSc=False,iplot=False,useController=True,htmlplot=False,outdir='./_output'):
    """
    Example python script for solving the 2d acoustics equations.
    """

    if use_PETSc:
        from petclaw.grid import Dimension, Grid
        from petclaw.evolve.clawpack import PetClawSolver2D as mySolver
        output_format='petsc'
    else:
        from pyclaw.grid import Dimension, Grid
        from pyclaw.evolve.clawpack import ClawSolver2D as mySolver
        output_format='ascii'

    from pyclaw.solution import Solution
    from pyclaw.controller import Controller
    from petclaw import plot

    # Initialize grid
    mx=100; my=100
    x = Dimension('x',-2.0,2.0,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',-2.0,2.0,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])

    grid.meqn = 1
    grid.mbc = 2
    tfinal = 1.0

    if use_PETSc:
        grid.init_q_petsc_structures()

    qinit(grid)
    initial_solution = Solution(grid)

    solver = mySolver()
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 2
    solver.mthlim = [1]*solver.mwaves

    claw = Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.output_format = output_format
    claw.tfinal = tfinal
    claw.solutions['n'] = initial_solution
    claw.solver = solver
    claw.nout = 10

    # Solve
    status = claw.run()

    if htmlplot:  plot.plotHTML(outdir=outdir,format=output_format)
    if iplot:     plot.plotInteractive(outdir=outdir,format=claw.output_format)


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        kpp(*args,**kwargs)
    else: kpp()
