#!/usr/bin/env python
# encoding: utf-8
IC='gauss_square'
if IC=='gauss_square':
    beta=200.; x0=0.3; mx=100
elif IC=='wavepacket':
    beta=100.; x0=0.5; mx=100

def fig_61_62_63(solver_order=2, limiters=0):
    """
    Compare several methods for advecting a Gaussian and square wave.

    The settings coded here are for Figure 6.1(a).
    For Figure 6.1(b), set solver.order=2.
    For Figure 6.2(a), set solver.order=2 and solver.limiters = pyclaw.limiters.tvd.minmod (1)
    For Figure 6.2(b), set solver.order=2 and solver.limiters = pyclaw.limiters.tvd.superbee (2)
    For Figure 6.2(c), set solver.order=2 and solver.limiters = pyclaw.limiters.tvd.MC (4)

    For Figure 6.3, set IC='wavepacket' and other options as appropriate.
    """
    import numpy as np
    from clawpack import pyclaw
    from clawpack import riemann

    solver = pyclaw.ClawSolver1D(riemann.advection_1D)

    solver.bc_lower[0] = 2
    solver.bc_upper[0] = 2
    solver.limiters = limiters
    solver.order = solver_order
    solver.cfl_desired = 0.8

    x = pyclaw.Dimension(0.0,1.0,mx,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain, num_eqn)
    state.problem_data['u'] = 1.

    xc = domain.grid.x.centers
    if IC=='gauss_square':
        state.q[0,:] = np.exp(-beta * (xc-x0)**2) + (xc>0.6)*(xc<0.8)
    elif IC=='wavepacket':
        state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.sin(80.*xc)
    else:
        raise Exception('Unrecognized initial condition specification.')

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver

    claw.tfinal = 10.0
    return claw


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(fig_61_62_63)
