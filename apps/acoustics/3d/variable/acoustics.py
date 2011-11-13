#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def acoustics3D(iplot=False,htmlplot=False,use_petsc=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 3d acoustics equations.
    """

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver3D()
    else:
        raise Exception('Unrecognized solver_type.')

    solver.dim_split=False
    solver.mwaves = 3
    solver.limiters = pyclaw.limiters.tvd.MC

    solver.bc_lower[0]=pyclaw.BC.reflecting
    solver.bc_upper[0]=pyclaw.BC.outflow
    solver.bc_lower[1]=pyclaw.BC.reflecting
    solver.bc_upper[1]=pyclaw.BC.outflow
    solver.bc_lower[2]=pyclaw.BC.reflecting
    solver.bc_upper[2]=pyclaw.BC.outflow

    solver.aux_bc_lower[0]=pyclaw.BC.reflecting
    solver.aux_bc_upper[0]=pyclaw.BC.outflow
    solver.aux_bc_lower[1]=pyclaw.BC.reflecting
    solver.aux_bc_upper[1]=pyclaw.BC.outflow
    solver.aux_bc_lower[2]=pyclaw.BC.reflecting
    solver.aux_bc_upper[2]=pyclaw.BC.outflow

    solver.dim_split=True

    # Initialize grid
    mx=30; my=30; mz=30
    x = pyclaw.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.Dimension('y',-1.0,1.0,my)
    z = pyclaw.Dimension('z',-1.0,1.0,mz)
    grid = pyclaw.Grid([x,y,z])

    meqn = 4
    maux = 2 # density, sound speed
    state = pyclaw.State(grid,meqn,maux)

    zl = 2.0  # Impedance in left half
    cl = 2.0  # Sound speed in left half
    zr = 1.0  # Impedance in right half
    cr = 1.0  # Sound speed in right half

    grid.compute_c_center()
    X,Y,Z = grid._c_center

    state.aux[0,:,:,:] = zl*(X<0.) + zr*(X>=0.) # Impedance
    state.aux[1,:,:,:] = cl*(X<0.) + cr*(X>=0.) # Sound speed

    x0 = -0.5; y0 = 0.; z0 = 0.
    r = np.sqrt((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)
    width=0.2
    state.q[0,:,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:,:] = 0.
    state.q[2,:,:,:] = 0.

    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir=outdir

    # Solve
    claw.tfinal = 0.6
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir,format=claw.output_format)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir,format=claw.output_format)

if __name__=="__main__":
    import sys
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics3D)
