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
    """
    from numpy import sin, cos
    from copy import deepcopy    

    # Polar coordinates (x coordinate = radius,  y coordinate = theta)
    nbrCells = len(mC[0])
    pC = deepcopy(mC)

    for iC in range(0,nbrCells):
        pC[0][iC] = mC[0][iC] * cos(mC[1][iC])
        pC[1][iC] = mC[0][iC] * sin(mC[1][iC])
    
    
    return pCCenters


def qinit(state):
    
    # Compute location of all grid cell center coordinates and store them
    state.grid.compute_p_center(recompute=True)

   
def setaux(state):
    """ 
    Set auxiliary array
    """    
    
    # Compute location of all grid cell edge coordinates and store them
    state.grid.compute_p_edge(recompute=True)

    # Define the auxiliary vector
    mx = len(state.grid.x.center) # Number of cell in the x direction
    my = len(state.grid.y.center) # Number of cell in the y direction 
    aux = np.empty((3,mx,my), order='F')

    # Define array that will contain the four nodes of a cell
    xp = np.zeros([5,1])
    yp = np.zeros([5,1])

    # Get grid spacing
    dxc = state.grid.d[0]
    dyc = state.grid.d[1]

    # Set auxiliary array
    # aux[1,i,j] is edge velocity at "left" boundary of grid point (i,j)
    # aux[2,i,j,] is edge velocity at "bottom" boundary of grid point (i,j)
    # aux[3,i,j] = kappa  is ratio of cell area to (dxc * dyc)
    for i in range(0,mx):
        for j in range(0,my):
            xp[0] = state.grid.p_edge[0][i][j]
            yp[0] = state.grid.p_edge[1][i][j]

            xp[1] = state.grid.p_edge[0][i][j+1]
            yp[1] = state.grid.p_edge[1][i][j+1]

            xp[2] = state.grid.p_edge[0][i+1][j+1]
            yp[2] = state.grid.p_edge[1][i+1][j+1]

            xp[3] = state.grid.p_edge[0][i+1][j]
            yp[3] = state.grid.p_edge[1][i+1][j]

            aux[0,i,j] = (stream(xp[1],yp[1])- stream(xp[0],yp[0]))/dyc
            aux[1,i,j] = -(stream(xp[3],yp[3])- stream(xp[0],yp[0]))/dxc

            xp[4] = xp[0]
            yp[4] = yp[0]
            area = 0.

            for iNode in range(0,4):
                area = area + 1./2.*(yp[iNode]+yp[iNode+1])*(xp[iNode+1]-xp[iNode])
        	    
            aux[2,i,j] = area/(dxc*dyc)
            


def stream(xp,yp):
    """ 
    Calculates the stream function in physical space.
    Clockwise rotation. One full rotation corresponds to 1 (second).
    """
    streamValue = np.pi*(np.square(xp) + np.square(yp))

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

    solver.mwaves = 1

    solver.dim_split = 0

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
    mx = 4

    ylower = 0.0
    yupper = np.pi*2.0
    my = 6

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])
    #grid.mapc2p = mapc2p_annulus # Override the default mapc2p function implemented in grid.py


    # Sate:
    state = pyclaw.State(grid)
    state.meqn = 1  # Number of equations

    # Set initial solution
    # ====================
    # qinit(state) # This function is defined above

    # Set auxiliary array
    # ===================
    state.aux = setaux(state)

  



if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        error=advection_annulus(*args,**kwargs)
        print 'Error: ',error
    else: advection_annulus()





