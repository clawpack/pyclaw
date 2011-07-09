#!/usr/bin/env python
# encoding: utf-8

#===========================================================================
# Import libraries
#===========================================================================
import numpy as np
from petsc4py import PETSc
import petclaw


def mapc2p_annulus(grid,mCCenters):
    """
    Specifies the mapping to curvilinear coordinates.

    Takes as input:    array_list made by x_coordinates, y_ccordinates in the map space
    Returns as output: array_list made by x_coordinates, y_ccordinates in the physical space
    """
    from numpy import sin, cos
    from copy import deepcopy    

    # Polar coordinates (x coordinate = radius,  y coordinate = theta)
    nbrCC = len(mCCenters[0])
    pCCenters = deepcopy(mCCenters)

    for iCC in range (0,nbrCC):
        pCCenters[0][iCC] = mCCenters[0][iCC] * cos(mCCenters[1][iCC])
        pCCenters[1][iCC] = mCCenters[0][iCC] * sin(mCCenters[1][iCC])
    
    
    return pCCenters


def qinit(state):
    
    x = state.grid.x.center
    y = state.grid.y.center
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] > 0.0 and x[i] < 0.5 and y[j]>0.0 and y[j] < 0.5:
                state.q[:,i,j] = 1.0
            else:
                state.q[:,i,j] = 0.1



def advection_annulus(use_petsc=False,iplot=0,htmlplot=False,outdir='./_output',solver_type='classic'):
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np

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
    # Initialize grids, then initialize the solution associated to the grid and
    # finally initialize aux array
    #===========================================================================

    # Grid:
    xlower = 0.2
    xupper = 1.0
    mx = 40

    ylower = 0.0
    yupper = np.pi*2.0
    my=120

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])



