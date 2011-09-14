#!/usr/bin/env python
# encoding: utf-8
"""
Example python script for solving the 1d variable-coefficient advection 
equation: q_t + u(x)q_x = 0.
"""

import numpy as np

def qinit(state):

    # Initial Data parameters
    ic = 3
    beta = 100.
    gamma = 0.
    x0 = 0.3
    x1 = 0.7
    x2 = 0.9

    x =state.grid.x.center
    
    # Gaussian
    qg = np.exp(-beta * (x-x0)**2) * np.cos(gamma * (x - x0))
    # Step Function
    qs = (x > x1) * 1.0 - (x > x2) * 1.0
    
    if   ic == 1: state.q[0,:] = qg
    elif ic == 2: state.q[0,:] = qs
    elif ic == 3: state.q[0,:] = qg + qs


def auxinit(state):
    # Initilize petsc Structures for aux
    xc=state.grid.x.center
    state.aux[0,:] = np.sin(2.*np.pi*xc)+2
    

def vc_advection(use_petsc=False,solver_type='classic',kernel_language='Python',iplot=False,htmlplot=False,outdir='./_output'):

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()

    solver.kernel_language = kernel_language
    from riemann import rp_vc_advection
    solver.mwaves = rp_vc_advection.mwaves
    if solver.kernel_language=='Python': 
        solver.rp = rp_vc_advection.rp_vc_advection_1d
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.mthbc_lower[0] = 2
    solver.mthbc_upper[0] = 2
    solver.mthauxbc_lower[0] = 2
    solver.mthauxbc_upper[0] = 2

    xlower=0.0; xupper=1.0; mx=100
    x    = pyclaw.Dimension('x',xlower,xupper,mx)
    grid = pyclaw.Grid(x)
    maux=1
    meqn = 1
    state = pyclaw.State(grid,meqn,maux)

    qinit(state)
    auxinit(state)

    claw = pyclaw.Controller()
    claw.outdir = outdir
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver

    claw.tfinal = 1.0
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(vc_advection)
