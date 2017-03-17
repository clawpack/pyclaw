#!/usr/bin/env python
# encoding: utf-8
from __future__ import absolute_import
IC='gauss_square'
if IC=='gauss_square':
    beta=200.; x0=0.3; mx=100
elif IC=='wavepacket':
    beta=100.; x0=0.5; mx=100

def fig_61_62_63(kernel_language='Python',iplot=False,htmlplot=False,solver_type='classic',outdir='./_output'):
    """
    Compare several methods for advecting a Gaussian and square wave.

    The settings coded here are for Figure 6.1(a).
    For Figure 6.1(b), set solver.order=2.
    For Figure 6.2(a), set solver.order=2 and solver.limiters = pyclaw.limiters.tvd.minmod
    For Figure 6.2(b), set solver.order=2 and solver.limiters = pyclaw.limiters.tvd.superbee
    For Figure 6.2(c), set solver.order=2 and solver.limiters = pyclaw.limiters.tvd.MC

    For Figure 6.3, set IC='wavepacket' and other options as appropriate.
    """
    import numpy as np
    from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()

    solver.kernel_language = kernel_language
    from clawpack.riemann import rp_advection
    solver.num_waves = rp_advection.num_waves
    if solver.kernel_language=='Python': 
        solver.rp = rp_advection.rp_advection_1d

    solver.bc_lower[0] = 2
    solver.bc_upper[0] = 2
    solver.limiters = 0
    solver.order = 1
    solver.cfl_desired = 0.8

    x = pyclaw.Dimension('x',0.0,1.0,mx)
    grid = pyclaw.Grid(x)
    num_eqn = 1
    state = pyclaw.State(grid,num_eqn)
    state.problem_data['u']=1.

    xc=grid.x.center
    if IC=='gauss_square':
        state.q[0,:] = np.exp(-beta * (xc-x0)**2) + (xc>0.6)*(xc<0.8)
    elif IC=='wavepacket':
        state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.sin(80.*xc)
    else:
        raise Exception('Unrecognized initial condition specification.')

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir

    claw.tfinal =10.0
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(fig_61_62_63)
