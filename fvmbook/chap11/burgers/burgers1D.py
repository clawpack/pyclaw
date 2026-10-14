#!/usr/bin/env python
# encoding: utf-8

def burgers():
    """
    Example from Chapter 11 of LeVeque, Figure 11.8.
    Shows decay of an initial wave packet to an N-wave with Burgers' equation.
    """
    import numpy as np

    from clawpack import pyclaw
    from clawpack import riemann

    solver = pyclaw.ClawSolver1D(riemann.burgers_1D)

    solver.limiters = pyclaw.limiters.tvd.MC
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    x = pyclaw.Dimension(-8.0,8.0,1000,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    xc = domain.grid.x.centers
    state.q[0,:] = (xc>-np.pi)*(xc<np.pi)*(2.*np.sin(3.*xc)+np.cos(2.*xc)+0.2)
    state.q[0,:] = state.q[0,:]*(np.cos(xc)+1.)
    state.problem_data['efix']=True

    claw = pyclaw.Controller()
    claw.tfinal = 6.0
    claw.num_output_times   = 30
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(burgers)
