#!/usr/bin/env python
# encoding: utf-8

def acoustics(iplot=False,htmlplot=False,outdir='./_output'):
    """
    This example solves the 1-dimensional acoustics equations in a homogeneous
    medium.
    """
    import numpy as np
    from clawpack import riemann
    from clawpack import pyclaw

    solver = pyclaw.ClawSolver1D(riemann.acoustics_1D)

    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall

    solver.cfl_desired = 1.0
    solver.cfl_max     = 1.0

    x = pyclaw.Dimension(0.0,1.0,200,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)

    # Set problem-specific variables
    rho = 1.0
    bulk = 1.0
    state.problem_data['rho']=rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']=np.sqrt(rho*bulk)
    state.problem_data['cc']=np.sqrt(bulk/rho)

    # Set the initial condition
    xc = domain.grid.x.centers
    state.q[0,:] = np.cos(2*np.pi*xc)
    state.q[1,:] = 0.
    
    # Set up the controller
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.num_output_times = 40
    claw.tfinal = 2.0

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)
