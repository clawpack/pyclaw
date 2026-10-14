#!/usr/bin/env python
# encoding: utf-8
r"""
2D shallow water: radial dam break
==================================

Solve the 2D shallow water equations:

.. math::
    h_t + (hu)_x + (hv)_y = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y = 0 \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y = 0.

The initial condition is a circular area with high depth surrounded by lower-depth water.
The top and right boundary conditions reflect, while the bottom and left boundaries
are outflow.
"""

import numpy as np
from clawpack import riemann
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn

def qinit(state,h_in=2.,h_out=1.,dam_radius=0.5):
    x0=0.
    y0=0.
    X, Y = state.p_centers
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)

    state.q[depth     ,:,:] = h_in*(r<=dam_radius) + h_out*(r>dam_radius)
    state.q[x_momentum,:,:] = 0.
    state.q[y_momentum,:,:] = 0.

    
def setup(kernel_language='Fortran', use_petsc=False, outdir='./_output',
          solver_type='classic', riemann_solver='roe',disable_output=False):
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if riemann_solver.lower() == 'roe':
        rs = riemann.shallow_roe_with_efix_2D
    elif riemann_solver.lower() == 'hlle':
        rs = riemann.shallow_hlle_2D

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(rs)
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=1
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(rs)

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.wall

    # Domain:
    xlower = -2.5
    xupper = 2.5
    mx = 150
    ylower = -2.5
    yupper = 2.5
    my = 150
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0

    qinit(state)

    claw = pyclaw.Controller()
    claw.tfinal = 2.5
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    if disable_output:
        claw.output_format = None
    claw.outdir = outdir
    claw.num_output_times = 10
    claw.setplot = setplot
    claw.keep_copy = True

    return claw

#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2.5, 2.5]
    plotaxes.ylimits = [-2.5, 2.5]
    plotaxes.title = 'Water height'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.5
    plotitem.pcolor_cmax = 1.5
    plotitem.add_colorbar = True
    
    # Scatter plot of depth
    plotfigure = plotdata.new_plotfigure(name='Scatter plot of h', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0., 2.5]
    plotaxes.ylimits = [0., 2.1]
    plotaxes.title = 'Scatter plot of h'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = depth
    def q_vs_radius(current_data):
        from numpy import sqrt
        x = current_data.x
        y = current_data.y
        r = sqrt(x**2 + y**2)
        q = current_data.q[depth,:,:]
        return r,q
    plotitem.map_2d_to_1d = q_vs_radius
    plotitem.plotstyle = 'o'


    # Figure for x-momentum
    plotfigure = plotdata.new_plotfigure(name='Momentum in x direction', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2.5, 2.5]
    plotaxes.ylimits = [-2.5, 2.5]
    plotaxes.title = 'Momentum in x direction'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = x_momentum
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = False       # show on plot?
    

    # Figure for y-momentum
    plotfigure = plotdata.new_plotfigure(name='Momentum in y direction', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2.5, 2.5]
    plotaxes.ylimits = [-2.5, 2.5]
    plotaxes.title = 'Momentum in y direction'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = y_momentum
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = False       # show on plot?
    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
