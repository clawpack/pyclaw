#!/usr/bin/env python
# encoding: utf-8
r"""
Two-dimensional advection
=========================

Solve the two-dimensional linear advection equation

.. math:: 
    q_t + (uq)_x + (vq)_y & = 0

Here q is a conserved quantity, and (u,v) is the velocity vector.
"""

import numpy as np

def qinit(state):

    # Set initial conditions for q.
    # Sample scalar equation with data that is piecewise constant with
    # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    #     0.1  otherwise
    
    x = state.grid.x.centers
    y = state.grid.y.centers
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] > 0.0 and x[i] < 0.5 and y[j]>0.0 and y[j] < 0.5:
                state.q[:,i,j] = 1.0
            else:
                state.q[:,i,j] = 0.1
                
def setup(use_petsc=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d advection equation.
    """
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D(riemann.advection_2D)
        solver.dimensional_split = 1
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.advection_2D)

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic

    solver.cfl_max=1.0
    solver.cfl_desired = 0.9

    #===========================================================================
    # Initialize domain, then initialize the solution associated to the domain and
    # finally initialize aux array
    #===========================================================================

    # Domain:
    mx=50; my=50
    x = pyclaw.Dimension('x',0.0,1.0,mx)
    y = pyclaw.Dimension('y',0.0,1.0,my)
    domain = pyclaw.Domain([x,y])

    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['u'] = 0.5 # Parameters (global auxiliary variables)
    state.problem_data['v'] = 1.0

    # Initial solution
    # ================
    qinit(state) # This function is defined above


    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
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
