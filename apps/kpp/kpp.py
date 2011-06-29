#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(state,rad=1.0):
    x = state.grid.x.center
    y = state.grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    state.q[0,:,:] = 0.25*np.pi + 3.25*np.pi*(r<=rad)


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

    solver.mwaves = 1
    solver.mthbc_lower[0]=pyclaw.BC.outflow
    solver.mthbc_upper[0]=pyclaw.BC.outflow
    solver.mthbc_lower[1]=pyclaw.BC.outflow
    solver.mthbc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=200; my=200
    x = pyclaw.Dimension('x',-2.0,2.0,mx)
    y = pyclaw.Dimension('y',-2.0,2.0,my)
    grid = pyclaw.Grid([x,y])
    state = pyclaw.State(grid)
    state.meqn = 1

    qinit(state)

    solver.dim_split = 1
    solver.cfl_max = 1.0
    solver.cfl_desired = 0.9
    solver.mwaves = 2
    solver.limiters = pyclaw.limiters.tvd.minmod

    claw = pyclaw.Controller()
    claw.tfinal = 1.0
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.nout = 10

    # Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        kpp(*args,**kwargs)
    else: kpp()
