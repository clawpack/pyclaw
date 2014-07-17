#!/usr/bin/env python
# encoding: utf-8
"""
Compressible Euler flow in cylindrical symmetry
===============================================

Solve the Euler equations of compressible fluid dynamics in 2D r-z coordinates:

.. math::
    \rho_t + (\rho u)_x + (\rho v)_y & = - \rho v / r \\
    (\rho u)_t + (\rho u^2 + p)_x + (\rho uv)_y & = -\rho u v / r \\
    (\rho v)_t + (\rho uv)_x + (\rho v^2 + p)_y & = - \rho v^2 / r \\
    E_t + (u (E + p) )_x + (v (E + p))_y & = - (E + p) v / r.

Here :math:`\rho` is the density, (u,v) is the velocity, and E is the total energy.
The radial coordinate is denoted by r.

The problem involves a planar shock wave impacting a spherical low-density bubble.
The problem is 3-dimensional but has been reduced to two dimensions using 
cylindrical symmetry.

This problem demonstrates:

    - how to incorporate source (non-hyperbolic) terms using both Classic and SharpClaw solvers
    - how to impose a custom boundary condition
    - how to use the auxiliary array for spatially-varying coefficients
"""

import numpy as np
from clawpack import riemann

gamma = 1.4 # Ratio of specific heats
gamma1 = gamma - 1.
x0=0.5; y0=0.; r0=0.2

def incoming_shock(state,dim,t,qbc,num_ghost):
    """
    Incoming shock at left boundary.
    """
    rho = 1.4;
    ux = 3.0;
    uy = 0.0;
    p = 1.0;
    T = 2*p/rho;

    for i in xrange(num_ghost):
        qbc[0,i,...] = rho
        qbc[1,i,...] = rho*ux
        qbc[2,i,...] = rho*uy
        qbc[3,i,...] = rho*(5.0/4 * T + (ux**2+uy**2)/2.)


def setup(use_petsc=False,solver_type='classic', outdir='_output', kernel_language='Fortran',
        disable_output=False, mx=320, my=80, tfinal=2.0, num_output_times = 10):
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    import riemannsolver

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemannsolver)
    else:
        solver = pyclaw.ClawSolver2D(riemannsolver)
        solver.dimensional_split = True

    solver.num_eqn = 4
    solver.num_waves = 4
    x = pyclaw.Dimension('x',0.0,3.0,mx)
    y = pyclaw.Dimension('y',0.0,1.0,my)
    domain = pyclaw.Domain([x,y])

    num_aux=1
    state = pyclaw.State(domain,solver.num_eqn,num_aux)
    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.dt_initial=0.005
    solver.user_bc_lower = incoming_shock

    solver.bc_lower[0]=pyclaw.BC.custom
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.wall
    solver.bc_upper[1]=pyclaw.BC.wall
    #Aux variable in ghost cells doesn't matter
    solver.aux_bc_lower[0]=pyclaw.BC.extrap
    solver.aux_bc_upper[0]=pyclaw.BC.extrap
    solver.aux_bc_lower[1]=pyclaw.BC.wall
    solver.aux_bc_upper[1]=pyclaw.BC.wall

    rho = 1.4;
    ux = 3.0;
    uy = 0.;
    p = 1.0;
    T = 2*p/rho;

    state.q[0,...] = rho;
    state.q[1,...] = rho*ux;
    state.q[2,...] = rho*uy;
    state.q[3,...] = rho*(5.0/4 * T + (ux**2+uy**2)/2.);

    x,y = domain.grid.p_centers
    state.aux[0,...]= 1 - (x > 0.6)*(y < 0.2)

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.tfinal = tfinal
    claw.num_output_times = num_output_times
    claw.outdir = outdir
    claw.setplot = setplot

    return claw

    
def setplot(plotdata):
    """ 
    Plot solution using VisClaw.
    """ 
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Pressure plot
    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Density'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = fill_step
    #plotitem = plotaxes.new_plotitem(plot_type='2d_schlieren')
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.pcolor_cmin = 1
    plotitem.pcolor_cmax=7.0
    plotitem.plot_var = 0
    plotitem.add_colorbar = True
    

    # Energy plot
    plotfigure = plotdata.new_plotfigure(name='Energy', figno=2)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Energy'
    plotaxes.scaled = True      # so aspect ratio is 1
    
    plotaxes.afteraxes = fill_step
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.pcolor_cmin = 2.
    plotitem.pcolor_cmax=18.0
    plotitem.plot_var = 3
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    
    return plotdata
   
def fill_step(currentdata):
    import matplotlib.pyplot as plt
    plt.hold(True);
    rectangle = plt.Rectangle((0.6,0.0),2.4,0.2,color="k",fill=True)
    plt.gca().add_patch(rectangle)
    plt.hold(False);


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
