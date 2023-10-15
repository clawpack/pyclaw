#!/usr/bin/env python
# encoding: utf-8
r"""
Traffic flow
============

Solve the Lighthill-Whitham-Richards (LWR) traffic flow model:

.. math::
    q_t + u (q(1-q))_x & = 0.

Here q is the density of cars, and u is a constant specifying the speed limit.
"""

from __future__ import absolute_import
from clawpack import riemann

def setup(use_petsc=0,outdir='./_output',solver_type='classic'):
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann.traffic_1D)
    else:
        solver = pyclaw.ClawSolver1D(riemann.traffic_1D)

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    x = pyclaw.Dimension(-1.0,1.0,500,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    grid = state.grid
    xc=grid.p_centers[0]

    state.q[0,:] = 0.75*(xc<0) + 0.1*(xc>0.) 

    state.problem_data['efix']=True
    state.problem_data['umax']=1.

    claw = pyclaw.Controller()
    claw.tfinal =2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True

    return claw


def setplot(plotdata):
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.ylimits = [-0.1, 1.1]
    plotaxes.title = 'q[0]'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
