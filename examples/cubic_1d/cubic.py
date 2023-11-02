#!/usr/bin/env python
# encoding: utf-8

r"""
Cubic equation
=========================

Solve the cubic conservation law:

.. math::
    q_t + (q^3)_x = 0.

This is a scalar nonlinear conservation law which is often used as a simple
example for problems with non-convex flux functions.
"""
import numpy as np
from clawpack import riemann

def setup(use_petsc=0, outdir='./_output', solver_type='classic', weno_order=5, N=1000):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    riemann_solver = riemann.cubic_1D

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann_solver)
        solver.weno_order = weno_order
    else:
        solver = pyclaw.ClawSolver1D(riemann_solver)
        solver.limiters = pyclaw.limiters.tvd.vanleer

    solver.cfl_max = 1.0
    solver.cfl_desired = 0.5

    solver.kernel_language = 'Fortran'

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    x = pyclaw.Dimension(-1.0, 3.0, N, name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain, num_eqn)

    xc = state.grid.x.centers
    qL = 4.0
    qR = -2.0
    state.q[0,:] = (xc < -0.5) * qL + (xc >= -0.5) * qR

    claw = pyclaw.Controller()
    claw.tfinal = 0.2
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True

    return claw


def setplot(plotdata):
    """
    Plot solution using VisClaw.
    """
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'

    return plotdata


if __name__=="__main__":
    # Run the example from the command line.
    # You can use the usual options, e.g. `iplot=1` to plot
    # the numerical solution.
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup, setplot)
