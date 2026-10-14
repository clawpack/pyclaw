#!/usr/bin/env python
# encoding: utf-8
r"""
2D shallow water: flow over a sill
==================================

Solve the 2D shallow water equations with
variable bathymetry:

.. :math:
    h_t + (hu)_x + (hv)_y & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y & = -g h b_x \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y & = -g h b_y.

The bathymetry contains a gaussian hump.
are outflow.
"""

from clawpack import riemann
from clawpack import pyclaw
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn
import numpy as np

amplitude = 0.1  # Height of incoming wave
t_bdy = 5.       # Stop sending in waves at this time

def bathymetry(x,y):
    r2 = (x-1.)**2 + (y-0.5)**2
    return 0.8*np.exp(-10*r2)

def wave_maker_bc(state,dim,t,qbc,auxbc,num_ghost):
    "Generate waves at left boundary as if there were a moving wall there."
    if dim.on_lower_boundary:
        qbc[0,:num_ghost,:]=qbc[0,num_ghost,:]
        t=state.t
        if t <= t_bdy:
            vwall = amplitude*(np.sin(t*np.pi/1.5))
        else:
            vwall=0.
        for ibc in range(num_ghost-1):
            qbc[1,num_ghost-ibc-1,:] = 2*vwall - qbc[1,num_ghost+ibc,:]


def setup(kernel_language='Fortran', solver_type='classic', use_petsc=False,
          outdir='./_output'):

    solver = pyclaw.ClawSolver2D(riemann.shallow_bathymetry_fwave_2D)
    solver.dimensional_split = 1  # No transverse solver available

    solver.bc_lower[0] = pyclaw.BC.custom
    solver.user_bc_lower = wave_maker_bc
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.extrap
    solver.aux_bc_upper[1] = pyclaw.BC.extrap

    my = 20
    mx = 4*my
    x = pyclaw.Dimension(0.,4.,mx,name='x')
    y = pyclaw.Dimension(0,1, my,name='y')
    domain = pyclaw.Domain([x,y])
    state = pyclaw.State(domain,num_eqn,num_aux=1)

    X, Y = state.p_centers
    state.aux[0,:,:] = bathymetry(X,Y)

    state.q[depth,:,:] = 1. - state.aux[0,:,:]
    state.q[x_momentum,:,:] = 0.
    state.q[y_momentum,:,:] = 0.

    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1.e-3
    state.problem_data['sea_level'] = 0.

    claw = pyclaw.Controller()
    claw.tfinal = 10
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 40
    claw.setplot = setplot
    claw.keep_copy = True

    return claw

def surface_height(current_data):
    h = current_data.q[0,:,:]
    b = bathymetry(current_data.x,current_data.y)
    return h+b

def setplot(plotdata):
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Water height'
    plotaxes.scaled = False

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surface_height
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.9
    plotitem.pcolor_cmax = 1.1
    plotitem.add_colorbar = True

    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
