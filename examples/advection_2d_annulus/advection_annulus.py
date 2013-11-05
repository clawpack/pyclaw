#!/usr/bin/env python
# encoding: utf-8
r"""
Advection in an annular domain
==============================

Solve the linear advection equation:

.. math:: 
    q_t + (u(x,y) q)_x + (v(x,y) q)_y & = 0.

Here q is the density of some conserved quantity and (u,v) is the velocity
field.  We take a rotational velocity field: :math:`u = \cos(\theta), v = \sin(\theta)`.

This is the simplest example that shows how to use a mapped grid in PyClaw.
"""


#===========================================================================
# Import libraries
#===========================================================================
import numpy as np

def mapc2p_annulus(grid,mC):
    """
    Specifies the mapping to curvilinear coordinates.

    Takes as input: array_list made by x_coordinates, y_ccordinates in the map 
                    space.
    Returns as output: array_list made by x_coordinates, y_ccordinates in the 
                       physical space.

    Inputs: mC = list composed by two arrays 
                 [array ([xc1, xc2, ...]), array([yc1, yc2, ...])]

    Output: pC = list composed by two arrays 
                 [array ([xp1, xp2, ...]), array([yp1, yp2, ...])]
    """  
    # Define new empty list
    pC = []

    # Populate it with the physical coordinates 
    # Polar coordinates (x coordinate = radius,  y coordinate = theta)
    pC.append(mC[0][:]*np.cos(mC[1][:]))
    pC.append(mC[0][:]*np.sin(mC[1][:]))
    
    return pC


def qinit(state,mx,my):
    """
    Initialize with two Gaussian pulses.
    """

    # The following parameters match the vaules used in clawpack
    # ==========================================================
    # First gaussian pulse
    A1    = 1.    # Amplitude
    beta1 = 40.   # Decay factor
    x1    = -0.5  # x-coordinate of the centers
    y1    = 0.    # y-coordinate of the centers

    # Second gaussian pulse
    A2    = -1.   # Amplitude
    beta2 = 40.   # Decay factor
    x2    = 0.5   # x-coordinate of the centers
    y2    = 0.    # y-coordinate of the centers

    
    # Compute location of all grid cell centers coordinates and store them
    state.grid.compute_p_centers(recompute=True)

    xp = state.grid.p_centers[0]
    yp = state.grid.p_centers[1]
    state.q[0,:,:] = A1*np.exp(-beta1*(np.square(xp-x1) + np.square(yp-y1)))\
                   + A2*np.exp(-beta2*(np.square(xp-x2) + np.square(yp-y2)))


def setaux(state,mx,my):
    """ 
    Set auxiliary array
    aux[0,i,j] is edges velocity at "left" boundary of grid point (i,j)
    aux[1,i,j] is edges velocity at "bottom" boundary of grid point (i,j)
    aux[2,i,j] = kappa  is ratio of cell area to (dxc * dyc)
    """    
    
    # Compute location of all grid cell corner coordinates and store them
    state.grid.compute_p_edges(recompute=True)

    # Get grid spacing
    dxc = state.grid.delta[0]
    dyc = state.grid.delta[1]
    pcorners = state.grid.p_edges

    aux = velocities_capa(pcorners[0],pcorners[1],dxc,dyc)
    return aux


def velocities_upper(state,dim,t,auxbc,num_ghost):
    """
    Set the velocities for the ghost cells outside the outer radius of the annulus.
    """
    from mapc2p import mapc2p

    grid=state.grid
    mx = grid.num_cells[0]
    my = grid.num_cells[1]
    dxc = grid.delta[0]
    dyc = grid.delta[1]

    if dim == grid.dimensions[0]:
        xc1d = grid.lower[0]+dxc*(np.arange(mx+num_ghost,mx+2*num_ghost+1)-num_ghost)
        yc1d = grid.lower[1]+dyc*(np.arange(my+2*num_ghost+1)-num_ghost)
        yc,xc = np.meshgrid(yc1d,xc1d)

        xp,yp = mapc2p(xc,yc)

        auxbc[:,-num_ghost:,:] = velocities_capa(xp,yp,dxc,dyc)

    else:
        raise Exception('Custum BC for this boundary is not appropriate!')


def velocities_lower(state,dim,t,auxbc,num_ghost):
    """
    Set the velocities for the ghost cells outside the inner radius of the annulus.
    """
    from mapc2p import mapc2p

    grid=state.grid
    my = grid.num_cells[1]
    dxc = grid.delta[0]
    dyc = grid.delta[1]

    if dim == grid.dimensions[0]:
        xc1d = grid.lower[0]+dxc*(np.arange(num_ghost+1)-num_ghost)
        yc1d = grid.lower[1]+dyc*(np.arange(my+2*num_ghost+1)-num_ghost)
        yc,xc = np.meshgrid(yc1d,xc1d)

        xp,yp = mapc2p(xc,yc)

        auxbc[:,0:num_ghost,:] = velocities_capa(xp,yp,dxc,dyc)

    else:
        raise Exception('Custum BC for this boundary is not appropriate!')


def velocities_capa(xp,yp,dx,dy):

    mx = xp.shape[0]-1
    my = xp.shape[1]-1
    aux = np.empty((3,mx,my), order='F')

    # Bottom-left corners
    xp0 = xp[:mx,:my]
    yp0 = yp[:mx,:my]

    # Top-left corners
    xp1 = xp[:mx,1:]
    yp1 = yp[:mx,1:]

    # Top-right corners
    xp2 = xp[1:,1:]
    yp2 = yp[1:,1:]

    # Top-left corners
    xp3 = xp[1:,:my]
    yp3 = yp[1:,:my]

    # Compute velocity component
    aux[0,:mx,:my] = (stream(xp1,yp1)- stream(xp0,yp0))/dy
    aux[1,:mx,:my] = -(stream(xp3,yp3)- stream(xp0,yp0))/dx

    # Compute area of the physical element
    area = 1./2.*( (yp0+yp1)*(xp1-xp0) +
                   (yp1+yp2)*(xp2-xp1) +
                   (yp2+yp3)*(xp3-xp2) +
                   (yp3+yp0)*(xp0-xp3) )
    
    # Compute capa 
    aux[2,:mx,:my] = area/(dx*dy)

    return aux

    
def stream(xp,yp):
    """ 
    Calculates the stream function in physical space.
    Clockwise rotation. One full rotation corresponds to 1 (second).
    """
    streamValue = np.pi*(xp**2 + yp**2)

    return streamValue


def setup(use_petsc=False,outdir='./_output',solver_type='classic'):
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(riemann.vc_advection_2D)
        solver.dimensional_split = 0
        solver.transverse_waves = 2
        solver.order = 2
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.vc_advection_2D)

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic

    solver.aux_bc_lower[0] = pyclaw.BC.custom
    solver.aux_bc_upper[0] = pyclaw.BC.custom
    solver.user_aux_bc_lower = velocities_lower
    solver.user_aux_bc_upper = velocities_upper
    solver.aux_bc_lower[1] = pyclaw.BC.periodic
    solver.aux_bc_upper[1] = pyclaw.BC.periodic

    solver.dt_initial = 0.1
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.2

    solver.limiters = pyclaw.limiters.tvd.vanleer

    #===========================================================================
    # Initialize domain and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    # Domain:
    xlower = 0.2
    xupper = 1.0
    mx = 40

    ylower = 0.0
    yupper = np.pi*2.0
    my = 120

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    domain = pyclaw.Domain([x,y])
    domain.grid.mapc2p = mapc2p_annulus # Override default_mapc2p function implemented in geometry.py

    # State:
    num_eqn = 1  # Number of equations
    state = pyclaw.State(domain,num_eqn)

    
    # Set initial solution
    # ====================
    qinit(state,mx,my) # This function is defined above

    # Set auxiliary array
    # ===================
    state.aux = setaux(state,mx,my) # This function is defined above
    state.index_capa = 2

    
    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.output_style = 1
    claw.num_output_times = 10
    claw.tfinal = 1.0
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
    from mapc2p import mapc2p
    import numpy as np
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'
    plotaxes.afteraxes = "pylab.axis('scaled')" 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = -1.
    plotitem.pcolor_cmax = 1.
    plotitem.add_colorbar = True
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p


    # Figure for contour plot
    plotfigure = plotdata.new_plotfigure(name='contour', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_levels = np.linspace(-0.9, 0.9, 10)
    plotitem.contour_colors = 'k'
    plotitem.patchedges_show = 1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    return plotdata

 
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
