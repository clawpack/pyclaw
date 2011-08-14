#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def acoustics2D(iplot=False,htmlplot=False,use_petsc=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d acoustics equations.
    """

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver2D()
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D()

    solver.dim_split=False
    solver.mwaves = 2
    solver.limiters = pyclaw.limiters.tvd.MC

    solver.mthbc_lower[0]=pyclaw.BC.reflecting
    solver.mthbc_upper[0]=pyclaw.BC.outflow
    solver.mthbc_lower[1]=pyclaw.BC.reflecting
    solver.mthbc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=100; my=100
    x = pyclaw.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.Dimension('y',-1.0,1.0,my)
    grid = pyclaw.Grid([x,y])
    state = pyclaw.State(grid)

    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    state.aux_global['rho']= rho
    state.aux_global['bulk']=bulk
    state.aux_global['zz']= zz
    state.aux_global['cc']=cc

    state.meqn = 3

    Y,X = np.meshgrid(grid.y.center,grid.x.center)
    r = np.sqrt(X**2 + Y**2)
    width=0.2
    state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir=outdir

    # Solve
    claw.tfinal = 0.6
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir,format=claw.output_format)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir,format=claw.output_format)

    if use_petsc:
        pressure=claw.frames[claw.nout].state.gqVec.getArray().reshape([grid.ng[0],grid.ng[1],state.meqn])[:,:,0]
    else:
        pressure=claw.frames[claw.nout].state.q[:,:,0]
    return pressure


if __name__=="__main__":
    import sys
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics2D)
