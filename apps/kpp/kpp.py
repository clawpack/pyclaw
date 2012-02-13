#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(state,rad=1.0):
    x = state.grid.x.centers
    y = state.grid.y.centers
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    state.q[0,:,:] = 0.25*np.pi + 3.25*np.pi*(r<=rad)


def kpp(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d KPP equations.
    """

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()
    else:
        solver = pyclaw.ClawSolver2D()

    import riemann
    solver.rp = riemann.rp2_kpp

    solver.num_waves = 1
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.extrap
    solver.bc_upper[1]=pyclaw.BC.extrap

    # Initialize domain
    mx=200; my=200
    x = pyclaw.Dimension('x',-2.0,2.0,mx)
    y = pyclaw.Dimension('y',-2.0,2.0,my)
    domain = pyclaw.Domain([x,y])
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    qinit(state)

    solver.dimensional_split = 1
    solver.cfl_max = 1.0
    solver.cfl_desired = 0.9
    solver.num_waves = 2
    solver.limiters = pyclaw.limiters.tvd.minmod

    claw = pyclaw.Controller()
    claw.tfinal = 1.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 10

    # Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(kpp)
