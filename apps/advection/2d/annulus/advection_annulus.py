#!/usr/bin/env python
# encoding: utf-8

#===========================================================================
# Import libraries
#===========================================================================
import numpy as np
from petsc4py import PETSc
import petclaw


def mapc2p_annulus(grid,mC):
    """
    Specifies the mapping to curvilinear coordinates.

    Takes as input:    array_list made by x_coordinates, y_ccordinates in the map space
    Returns as output: array_list made by x_coordinates, y_ccordinates in the physical space

    Inputs: mC = list composed by two array [array ([xc1, xc2, ...]), array([yc1, yc2, ...])]

    Output: pC = list composed by two array [array ([xp1, xp2, ...]), array([yp1, yp2, ...])]
    """  

    # Polar coordinates (x coordinate = radius,  y coordinate = theta)
    nbrCells = len(mC[0])

    # Define new empty list
    pC = []

    # Populate it with the physical coordinates 
    pC.append(mC[0][:]*np.cos(mC[1][:]))
    pC.append(mC[0][:]*np.sin(mC[1][:]))
    
    return pC


def qinit(state,mx,my):
    """
    Computes the initial condition.
    """

    # The following parameters match the vaules used in clawpack
    # ==========================================================
    # First gaussian pulse
    A1    = 1.    # Amplitude
    beta1 = 40.   # Decay factor
    x1    = -0.5  # x-coordinate of the center
    y1    = 0.    # y-coordinate of the center

    # Second gaussian pulse
    A2    = -1.   # Amplitude
    beta2 = 40.   # Decay factor
    x2    = 0.5   # x-coordinate of the center
    y2    = 0.    # y-coordinate of the center

    
    # Compute location of all grid cell center coordinates and store them
    state.grid.compute_p_center(recompute=True)

    xp = state.grid.p_center[0]
    yp = state.grid.p_center[1]
    state.q[0,:,:] = A1*np.exp(-beta1*(np.square(xp-x1) + np.square(yp-y1)))\
                   + A2*np.exp(-beta2*(np.square(xp-x2) + np.square(yp-y2)))


    #print state.q[0,0,0],state.q[0,0,1],state.q[0,0,2],state.q[0,0,3],state.q[0,0,4]
    #print state.q[0,1,0],state.q[0,1,1],state.q[0,1,2],state.q[0,1,3],state.q[0,1,4]
    #print state.q[0,2,0],state.q[0,2,1],state.q[0,2,2],state.q[0,2,3],state.q[0,2,4]
    #print state.q[0,3,0],state.q[0,3,1],state.q[0,3,2],state.q[0,3,3],state.q[0,3,4]
    #print state.q[0,4,0],state.q[0,4,1],state.q[0,4,2],state.q[0,4,3],state.q[0,4,4]

   
def setaux(state,mx,my):
    """ 
    Set auxiliary array
    """    
    
    # Compute location of all grid cell edge coordinates and store them
    state.grid.compute_p_edge(recompute=True)

    # Define the auxiliary array
    aux = np.empty((3,mx,my), order='F')

    # Define array that will contain the four nodes of a cell
    xp = np.zeros([5,1])
    yp = np.zeros([5,1])

    # Get grid spacing
    dxc = state.grid.d[0]
    dyc = state.grid.d[1]

    # Set auxiliary array
    # aux[0,i,j] is edge velocity at "left" boundary of grid point (i,j)
    # aux[1,i,j] is edge velocity at "bottom" boundary of grid point (i,j)
    # aux[2,i,j] = kappa  is ratio of cell area to (dxc * dyc)
    for j in range(my):
        for i in range(mx):
            xp[0] = state.grid.p_edge[0][i,j]
            yp[0] = state.grid.p_edge[1][i,j]

            xp[1] = state.grid.p_edge[0][i,j+1]
            yp[1] = state.grid.p_edge[1][i,j+1]

            xp[2] = state.grid.p_edge[0][i+1,j+1]
            yp[2] = state.grid.p_edge[1][i+1,j+1]

            xp[3] = state.grid.p_edge[0][i+1,j]
            yp[3] = state.grid.p_edge[1][i+1,j]

            aux[0,i,j] = (stream(xp[1],yp[1])- stream(xp[0],yp[0]))/dyc
            aux[1,i,j] = -(stream(xp[3],yp[3])- stream(xp[0],yp[0]))/dxc

            xp[4] = xp[0]
            yp[4] = yp[0]
            area = 0.

            for iNode in range(4):
                area = area + 1./2.*(yp[iNode]+yp[iNode+1])*(xp[iNode+1]-xp[iNode]) 
            
            aux[2,i,j] = area/(dxc*dyc)

    #import matplotlib.pyplot as plt
    #plt.quiver(state.grid.p_center[0],state.grid.p_center[1],aux[0,:,:],aux[1,:,:])
    #plt.draw()

    #print aux[1,0,0],aux[1,0,1],aux[1,0,2],aux[1,0,3],aux[1,0,4]
    #print aux[1,1,0],aux[1,1,1],aux[1,1,2],aux[1,1,3],aux[1,1,4]
    #print aux[1,2,0],aux[1,2,1],aux[1,2,2],aux[1,2,3],aux[1,2,4]
    #print aux[1,3,0],aux[1,3,1],aux[1,3,2],aux[1,3,3],aux[1,3,4]
    #print aux[1,4,0],aux[1,4,1],aux[1,4,2],aux[1,4,3],aux[1,4,4]

    return aux
  

def stream(xp,yp):
    """ 
    Calculates the stream function in physical space.
    Clockwise rotation. One full rotation corresponds to 1 (second).
    """
    streamValue = np.pi*(xp**2 + yp**2)

    return streamValue


def advection_annulus(use_petsc=False,iplot=0,htmlplot=False,outdir='./_output',solver_type='classic'):
    #===========================================================================
    # Import libraries
    #===========================================================================
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D()
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D()


    solver.mthbc_lower[0] = pyclaw.BC.outflow
    solver.mthbc_upper[0] = pyclaw.BC.outflow
    solver.mthbc_lower[1] = pyclaw.BC.periodic
    solver.mthbc_upper[1] = pyclaw.BC.periodic

    solver.mthauxbc_lower[0] = pyclaw.BC.outflow
    solver.mthauxbc_upper[0] = pyclaw.BC.outflow
    solver.mthauxbc_lower[1] = pyclaw.BC.periodic
    solver.mthauxbc_upper[1] = pyclaw.BC.periodic

    solver.mwaves = 1

    solver.dim_split = 0
    solver.order_trans = 2
    solver.order = 2

    solver.dt_initial = 0.1
    solver.cfl_max = 1.0
    solver.cfl_desired = 0.9

    solver.limiters = pyclaw.limiters.tvd.vanleer

    #===========================================================================
    # Initialize grid and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    # Grid:
    xlower = 0.2
    xupper = 1.0
    mx = 40

    ylower = 0.0
    yupper = np.pi*2.0
    my = 120

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])
    grid.mapc2p = mapc2p_annulus # Override default_mapc2p function implemented in grid.py

    #print grid.d[0],grid.d[1]


    # Sate:
    state = pyclaw.State(grid)
    state.meqn = 1  # Number of equations

    
    # Set initial solution
    # ====================
    qinit(state,mx,my) # This function is defined above

    # Set auxiliary array
    # ===================
    state.aux = setaux(state,mx,my) # This function is defined above
    state.mcapa = 2

    
    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.outstyle = 1
    claw.nout = 20
    claw.tfinal = 0.5
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir

    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    #===========================================================================
    # Plot results
    #===========================================================================
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)



if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        error=advection_annulus(*args,**kwargs)
        print 'Error: ',error
    else: advection_annulus()





