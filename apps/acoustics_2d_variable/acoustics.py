#!/usr/bin/env python
# encoding: utf-8
r"""
Variable-coefficient acoustics example.
"""
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

    solver.bc_lower[0]=pyclaw.BC.reflecting
    solver.bc_upper[0]=pyclaw.BC.outflow
    solver.bc_lower[1]=pyclaw.BC.reflecting
    solver.bc_upper[1]=pyclaw.BC.outflow
    solver.aux_bc_lower[0]=pyclaw.BC.reflecting
    solver.aux_bc_upper[0]=pyclaw.BC.outflow
    solver.aux_bc_lower[1]=pyclaw.BC.reflecting
    solver.aux_bc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=200; my=200
    x = pyclaw.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.Dimension('y',-1.0,1.0,my)
    grid = pyclaw.Grid([x,y])

    meqn = 3
    maux = 2 # density, sound speed
    state = pyclaw.State(grid,meqn,maux)

    # Cell center coordinates
    Y,X = np.meshgrid(grid.y.center,grid.x.center)

    # Set aux arrays
    rhol = 4.0
    rhor = 1.0
    bulkl = 4.0
    bulkr = 4.0
    cl = np.sqrt(bulkl/rhol)
    cr = np.sqrt(bulkr/rhor)
    state.aux[0,:,:] = rhol*(X<0.) + rhor*(X>=0.) # Density
    state.aux[1,:,:] = cl*(X<0.) + cr*(X>=0.) # Sound speed

    # Set initial condition
    x0 = -0.5; y0 = 0.
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    width=0.1; rad=0.25
    state.q[0,:,:] = (np.abs(r-rad)<=width)*(1.+np.cos(np.pi*(r-rad)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir=outdir
    claw.nout = 20

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
