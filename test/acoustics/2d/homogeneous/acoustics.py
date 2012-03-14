#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(state,width=0.2):
    
    grid = state.grid
    x =grid.x.centers
    y =grid.y.centers
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.


def acoustics2D(use_petsc=False,kernel_language='Fortran',iplot=False,htmlplot=False,solver_type='classic', outdir = './_output', num_output_times = 10):
    """
    Example python script for solving the 2d acoustics equations.
    """
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D()
        solver.dimensional_split = 1
        solver.num_waves = 2
        solver.limiters = [4]*solver.num_waves
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()
        solver.num_waves = 2

    import riemann
    solver.rp = riemann.rp2_acoustics

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.extrap

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

    tfinal = 0.12

    qinit(state)
    initial_solution = pyclaw.Solution(state,domain)

    solver.dt_initial=np.min(domain.grid.delta)/state.problem_data['cc']*solver.cfl_desired

    claw = pyclaw.Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.tfinal = tfinal
    claw.solution = initial_solution
    claw.solver = solver
    claw.outdir = outdir
    claw.num_output_times = num_output_times

    # Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot()
    if iplot:     pyclaw.plot.interactive_plot()

    if use_petsc:
        grid = domain.grid
        pressure=claw.frames[claw.num_output_times].state.gqVec.getArray().reshape([state.num_eqn,grid.num_cells[0],grid.num_cells[1]],order='F')[0,:,:]
    else:
        pressure=claw.frames[claw.num_output_times].state.q[0,:,:]
    return pressure


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        acoustics2D(*args,**kwargs)
    else: acoustics2D()
