#!/usr/bin/env python
# encoding: utf-8
r"""
Two-dimensional acoustics
=========================

Solve the (linear) acoustics equations:

.. math:: 
    p_t + K (u_x + v_y) & = 0 \\ 
    u_t + p_x / \rho & = 0 \\
    v_t + p_y / \rho & = 0.

Here p is the pressure, (u,v) is the velocity, K is the bulk modulus,
and :math:`\rho` is the density.
"""
 
import numpy as np

def setup(kernel_language='Fortran',use_petsc=False,outdir='./_output',solver_type='classic',  disable_output=False):
    """
    Example python script for solving the 2d acoustics equations.
    """
    from clawpack import riemann
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver2D(riemann.acoustics_2D)
        solver.dimensional_split=True
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D(riemann.acoustics_2D)

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45

    solver.limiters = pyclaw.limiters.tvd.MC

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.extrap
    solver.bc_upper[1]=pyclaw.BC.extrap

    # Initialize domain
    mx=100; my=100
    x = pyclaw.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.Dimension('y',-1.0,1.0,my)
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    state = pyclaw.State(domain,num_eqn)

    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    state.problem_data['rho']= rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']= zz
    state.problem_data['cc']=cc

    qinit(state)

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.solution = pyclaw.Solution(state,domain)
    solver.dt_initial=np.min(domain.grid.delta)/state.problem_data['cc']*solver.cfl_desired

    claw.solver = solver
    claw.outdir = outdir

    num_output_times = 10
    
    claw.num_output_times = num_output_times

    claw.tfinal = 0.12

    claw.setplot = setplot

    return claw

def qinit(state,width=0.2):
    
    grid = state.grid
    x =grid.x.centers
    y =grid.y.centers
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.



#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 

    import os
    if os.path.exists('./1drad/_output'):
        qref_dir = os.path.abspath('./1drad/_output')
    else:
        qref_dir = None
        print "Directory ./1drad/_output not found"

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True

    # Figure for scatter plot
    plotfigure = plotdata.new_plotfigure(name='scatter', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Scatter plot'

    # Set up for item on these axes: scatter of 2d data
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    
    def p_vs_r(current_data):
        # Return radius of each patch cell and p value in the cell
        from pylab import sqrt
        x = current_data.x
        y = current_data.y
        r = sqrt(x**2 + y**2)
        q = current_data.q
        p = q[0,:,:]
        return r,p

    plotitem.map_2d_to_1d = p_vs_r
    plotitem.plot_var = 0
    plotitem.plotstyle = 'ob'
    
    return plotdata

    
if __name__=="__main__":
    import sys
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
