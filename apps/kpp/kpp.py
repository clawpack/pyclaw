#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,rad=1.0):
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    grid.zeros_q()
    grid.q[0,:,:] = 0.25*np.pi + 3.25*np.pi*(r<=rad)


def kpp(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d acoustics equations.
    """

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()
    else:
        solver = pyclaw.ClawSolver2D()

    solver.mwaves = 5
    solver.mthbc_lower[0]=pyclaw.BC.outflow
    solver.mthbc_upper[0]=pyclaw.BC.outflow
    solver.mthbc_lower[1]=pyclaw.BC.outflow
    solver.mthbc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=100; my=50
    print mx, my
    x = pyclaw.Dimension('x',-2.0,2.0,mx)
    y = pyclaw.Dimension('y',-2.0,2.0,my)
    grid = pyclaw.Grid([x,y])
    grid.meqn = 1
    grid.mbc = solver.mbc

    qinit(grid)

    solver.dim_split = 1
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 2
    solver.mthlim = pyclaw.limiters.minmod

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = 1.0
    claw.solution = pyclaw.Solution(grid)
    claw.solver = solver
    claw.nout = 10

    # Solve
    status = claw.run()

    from pyclaw import plot
    if htmlplot:  plot.plotHTML(outdir=outdir)
    if iplot:     plot.plotInteractive(outdir=outdir)


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        kpp(*args,**kwargs)
    else: kpp()
