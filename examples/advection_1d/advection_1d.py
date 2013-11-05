#!/usr/bin/env python
# encoding: utf-8

r"""
One-dimensional advection
=========================

Solve the linear advection equation:

.. math:: 
    q_t + u q_x & = 0.

Here q is the density of some conserved quantity and u is the velocity.

The initial condition is a Gaussian and the boundary conditions are periodic.
The final solution is identical to the initial data because the wave has
crossed the domain exactly once.
"""

def setup(nx=100,kernel_language='Python', use_petsc=False,
          solver_type='classic', weno_order=5, outdir='./_output'):
    import numpy as np
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        if kernel_language == 'Fortran':
            solver = pyclaw.ClawSolver1D(riemann.advection_1D)
        elif kernel_language=='Python': 
            solver = pyclaw.ClawSolver1D(riemann.advection_1D_py.advection_1D)
    elif solver_type=='sharpclaw':
        if kernel_language == 'Fortran':
            solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
        elif kernel_language=='Python': 
            solver = pyclaw.SharpClawSolver1D(riemann.advection_1D_py.advection_1D)
        solver.weno_order=weno_order
    else: raise Exception('Unrecognized value of solver_type.')

    solver.kernel_language = kernel_language

    solver.bc_lower[0] = 2
    solver.bc_upper[0] = 2

    x = pyclaw.Dimension('x',0.0,1.0,nx)
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)
    state.problem_data['u']=1.

    grid = state.grid
    xc=grid.x.centers
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    if outdir is not None:
        claw.outdir = outdir
    else:
        claw.output_format = None

    claw.tfinal =1.0
    claw.setplot = setplot

    return claw

#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotfigure = plotdata.new_plotfigure(name='q', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.ylimits = [-.2,1.0]
    plotaxes.title = 'q'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':2,'markersize':5}
    
    return plotdata

 
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
