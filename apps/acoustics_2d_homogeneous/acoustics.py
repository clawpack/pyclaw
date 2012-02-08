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
        solver.dimensional_split=False
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D()

    from riemann import rp2_acoustics
    solver.rp = rp2_acoustics

    solver.num_waves = 2
    solver.limiters = pyclaw.limiters.tvd.MC

    solver.bc_lower[0]=pyclaw.BC.wall
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.wall
    solver.bc_upper[1]=pyclaw.BC.extrap

    # Initialize domain
    mx=100; my=100
    x = pyclaw.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.Dimension('y',-1.0,1.0,my)
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    state = pyclaw.State(domain,num_eqn)

    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    state.problem_data['rho']= rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']= zz
    state.problem_data['cc']=cc

    grid = state.grid
    Y,X = np.meshgrid(grid.y.centers,grid.x.centers)
    r = np.sqrt(X**2 + Y**2)
    width=0.2
    state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir=outdir

    # Solve
    claw.tfinal = 0.6
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir,file_format=claw.output_format)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir,file_format=claw.output_format)

    if use_petsc:
        pressure=claw.frames[claw.num_output_times].state.gqVec.getArray().reshape([grid.num_cells[0],grid.num_cells[1],state.num_eqn])[:,:,0]
    else:
        pressure=claw.frames[claw.num_output_times].state.q[:,:,0]
    return pressure


if __name__=="__main__":
    import sys
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics2D)
