#!/usr/bin/env python
# encoding: utf-8

r"""
One-dimensional advection
=========================

Solve the linear advection equation on a nonuniform grid:

.. math::
    q_t + u q_x = 0.

Here q is the density of some conserved quantity and u is the velocity.
Here we have a nonuniform grid, given by the transformation x**2 on grid [-0.5,0.5] to [-0.25,0.25].
The initial condition is a Gaussian centered at 0 and the boundary conditions are periodic.
The final solution is identical to the initial data because the wave has
crossed the domain exactly once, which takes computational time 0.5, because the speed is 1 and grid length 0.5.
"""

import numpy as np
from clawpack import riemann

def mapc2p_nonunif(xc):
    """This function takes the interval [-xL,xR] and squares the computational coordinate
       while keeping the negative coordinates to be squared and yet retain their negativity
       to the physical coordinate [-xL**2, xR**2]
    """
    neg = -1*(xc < 0) + (xc > 0)
    xp = xc**2
    xp = neg*xp
    return xp


def setup(nx=100, kernel_language='Fortran', use_petsc=False, solver_type='classic', weno_order=5,
          time_integrator='SSP104', outdir='./_output'):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    from clawpack import pyclaw
    import clawpack.pyclaw.geometry

    if kernel_language == 'Fortran':
        riemann_solver = riemann.advection_1D
    elif kernel_language == 'Python':
        riemann_solver = riemann.advection_1D_py.advection_1D

    # sharpclaw does not accommodate nonuniform grids
    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D(riemann_solver)
    else: raise Exception('Unrecognized value of solver_type.')

    solver.kernel_language = kernel_language
    solver.order = 1 
    solver.limiters = pyclaw.tvd.minmod
    solver.num_eqn=1
    solver.num_waves=1
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.aux_bc_lower[0] = pyclaw.BC.periodic
    solver.aux_bc_upper[0] = pyclaw.BC.periodic

    x = pyclaw.Dimension(-0.5,0.5,nx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,1,num_aux=1)
    state.problem_data['u'] = 1.  # Advection velocity
    state.index_capa = 0

    xc = state.grid.x.centers
    grid1d = state.grid

    # mapping to nonunif grid
    grid1d.mapc2p = mapc2p_nonunif
    state.aux = np.zeros((1,nx)) # capacity array dx_p/dx_c
    state.aux[0,:] = np.diff(grid1d.p_nodes)/np.diff(state.grid.x.nodes)

    # Initial data
    beta = 100; gamma = 0; x0 = 0.0
    state.q[0,:] = np.exp(-beta * (grid1d.p_centers-x0)**2) * np.cos(gamma * (grid1d.p_centers - x0))

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    claw.tfinal = 0.5 # one cycle
    claw.outdir = outdir
    claw.num_output_times = 10
    claw.nstepout = 1
    if outdir is None:
        claw.output_format = None

    claw.setplot = setplot

    return claw

def setplot(plotdata):
    """
    Plot solution using VisClaw.
    """
    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.mapc2p = mapc2p_nonunif
    plotfigure = plotdata.new_plotfigure(name='q', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-0.25,0.25]
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
