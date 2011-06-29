#!/usr/bin/env python
# encoding: utf-8
def advection(kernel_language='Python',iplot=False,htmlplot=False,use_petsc=False,solver_type='classic',outdir='./_output'):
    """
    Example python script for solving the 1d advection equation.
    """
    import numpy as np

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()

    solver.kernel_language = kernel_language
    from riemann import rp_advection
    solver.mwaves = rp_advection.mwaves
    if solver.kernel_language=='Python': 
        solver.rp = rp_advection.rp_advection_1d

    solver.mthbc_lower[0] = 2
    solver.mthbc_upper[0] = 2

    x = pyclaw.Dimension('x',0.0,1.0,100)
    grid = pyclaw.Grid(x)
    state = pyclaw.State(grid)
    state.aux_global['u']=1.
    state.meqn = rp_advection.meqn

    xc=grid.x.center
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir

    claw.tfinal =1.0
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    import sys
    from pyclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    output=advection(*args,**kwargs)
