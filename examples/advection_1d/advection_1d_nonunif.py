#!/usr/bin/env python
# encoding: utf-8

r"""
One-dimensional advection
=========================

Solve the linear advection equation:

.. math::
    q_t + u q_x = 0.

Here q is the density of some conserved quantity and u is the velocity.
The grid is nonuniform as shown by the aux variable, with first half of the grid being 0.9*dx and second half being 1.1*dx, 
where dx is the computational grid size i.e. 1/nx
The initial condition is a Gaussian and the boundary conditions are periodic.
The final solution is identical to the initial data because the wave has
crossed the domain exactly once.
"""
from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.pyclaw.plot import plot




def setup(nx=100, kernel_language='Python', use_petsc=False, solver_type='classic', weno_order=5,
          time_integrator='SSP104', outdir='./_output'):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Fortran':
        riemann_solver = riemann.advection_1D
    elif kernel_language == 'Python':
        riemann_solver = riemann.advection_1D_py.advection_1D

    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D(riemann_solver)
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann_solver)
        solver.weno_order = weno_order
        solver.time_integrator = time_integrator
        if time_integrator == 'SSPLMMk3':
            solver.lmm_steps = 5
            solver.check_lmm_cond = True
    else: raise Exception('Unrecognized value of solver_type.')

    solver.kernel_language = kernel_language
    solver.order=1
    solver.limiters = None
    solver.num_eqn=1
    solver.num_waves=1
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.aux_bc_lower[0] = pyclaw.BC.periodic
    solver.aux_bc_upper[0] = pyclaw.BC.periodic

    x = pyclaw.Dimension(0.0,1.0,nx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,1, num_aux=1)

    state.problem_data['u'] = 1.  # Advection velocity

    state.index_capa = 0
    state.aux = np.zeros((1,nx))
    state.aux[0, :int(nx)/2 + 1] = 0.9 
    state.aux[0, int(nx)/2 + 1:] = 1.1

    # Initial data
    xc = state.grid.x.centers
    delta_xc = (1.0-0.0)/nx
    xc[:int(nx/2)+1] -= 0.05*delta_xc
    xc[int(nx/2)+1:] += 0.05*delta_xc
 
    beta = 100; gamma = 0; x0 = 0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    claw.tfinal = 1.0
    claw.outdir = outdir
    if outdir is None:
        claw.output_format = None

    claw.run()


    claw.setplot = setplot

#    plot(setplot=setplot,outdir='./_output_hb',plotdir='./plots_hb',iplot=False, htmlplot=True)

  # return claw

def setplot(plotdata):
    """
    Plot solution using VisClaw.
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
    #from clawpack.pyclaw.util import run_app_from_main
    #output = run_app_from_main(setup,setplot)
    setup()
