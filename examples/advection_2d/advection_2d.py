#!/usr/bin/env python
# encoding: utf-8
r"""
Two-dimensional advection
=========================

Solve the two-dimensional advection equation

.. math::
    q_t + u q_x + v q_y = 0

Here q is a conserved quantity, and (u,v) is the velocity vector.
"""

import numpy as np
from clawpack import riemann


def qinit(state):
    """Set initial condition for q.
       Sample scalar equation with data that is piecewise constant with
       q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
           0.1  otherwise
    """
    X, Y = state.grid.p_centers
    state.q[0,:,:] = 0.9*(0.1<X)*(X<0.6)*(0.1<Y)*(Y<0.6) + 0.1
                

def setup(use_petsc=False,outdir='./_output',solver_type='classic',
          dimensional_split=1,transverse_waves=0):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(riemann.advection_2D)
        solver.dimensional_split = dimensional_split
        solver.transverse_waves = transverse_waves
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.advection_2D)

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic

    if not dimensional_split and not transverse_waves and solver_type == 'classic':
        solver.cfl_max = 0.5
        solver.cfl_desired = 0.45
    else:
        solver.cfl_max = 1.0
        solver.cfl_desired = 0.9

    # Domain:
    mx = 50; my = 50
    x = pyclaw.Dimension(0.0,1.0,mx,name='x')
    y = pyclaw.Dimension(0.0,1.0,my,name='y')
    domain = pyclaw.Domain([x,y])

    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['u'] = 0.5 # Advection velocity
    state.problem_data['v'] = 1.0

    qinit(state)

    claw = pyclaw.Controller()
    claw.tfinal = 2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    if outdir is None:
        claw.output_format = None
    claw.setplot = setplot
    claw.keep_copy = True

    return claw


def setplot(plotdata):
    """ 
    Plot solution using VisClaw.
    """ 
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar = True
    
    # Figure for contour plot
    plotfigure = plotdata.new_plotfigure(name='contour', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_nlevels = 20
    plotitem.contour_min = 0.01
    plotitem.contour_max = 0.99
    plotitem.amr_contour_colors = ['b','k','r']
    
    return plotdata

    
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
