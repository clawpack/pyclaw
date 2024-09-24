#!/usr/bin/env python
# encoding: utf-8

def fig_31_38(iplot=False,htmlplot=False,outdir='./_output'):
    r"""Produces the output shown in Figures 3.1 and 3.8 of the FVM book.
    These involve simple waves in the acoustics system."""
    from clawpack import pyclaw
    from clawpack import riemann
    import numpy as np

    solver = pyclaw.ClawSolver1D(riemann.acoustics_1D)

    solver.limiters = pyclaw.limiters.tvd.MC
    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.extrap

    x = pyclaw.Dimension(-1.0,1.0,800,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain, num_eqn)

    # Set problem-specific variables
    rho = 1.0
    bulk = 0.25
    state.problem_data['rho']=rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']=np.sqrt(rho*bulk)
    state.problem_data['cc']=np.sqrt(bulk/rho)

    # Set the initial condition
    xc = domain.grid.x.centers
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = 0.5*np.exp(-80 * xc**2) + 0.5*(np.abs(xc+0.2)<0.1)
    state.q[1,:] = 0.
    
    # Set up the controller
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = 3.0
    claw.num_output_times   = 30

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(fig_31_38)
