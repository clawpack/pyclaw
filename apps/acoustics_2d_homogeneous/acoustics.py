#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def acoustics2D(iplot=False,kernel_language='Fortran',htmlplot=False,use_petsc=False,outdir='./_output',solver_type='classic',  disable_output=False):
    """
    Example python script for solving the 2d acoustics equations.
    """

    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver2D()
        solver.dimensional_split=True
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D()

    if kernel_language != 'Fortran':
        raise Exception('Unrecognized value of kernel_language for 2D acoustics')

    from clawpack.riemann import rp2_acoustics
    solver.rp = rp2_acoustics

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45

    solver.num_waves = 2
    solver.limiters = pyclaw.limiters.tvd.MC

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.extrap
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

    qinit(state)

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.solution = pyclaw.Solution(state,domain)
    solver.dt_initial=np.min(domain.grid.delta)/state.problem_data['cc']*solver.cfl_desired

    claw.solver = solver
    claw.outdir = outdir

    num_output_times = 10
    
    claw.num_output_times = num_output_times

    # Solve
    claw.tfinal = 0.12
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir,file_format=claw.output_format)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir,file_format=claw.output_format)

    return claw.frames[-1].state

def qinit(state,width=0.2):
    
    grid = state.grid
    x =grid.x.centers
    y =grid.y.centers
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

if __name__=="__main__":
    import sys
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics2D)
