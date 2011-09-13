#!/usr/bin/env python
# encoding: utf-8


"""
2D shallow water equations on a spherical surface.
The approximation of the three-dimensional equations
is restricted to the surface of the sphere.
"""

import numpy as np
from petclaw import plot
#import pdb  # Debugger

Rsphere = 1.0
g = 1.0



def mapc2p_sphere(grid,mC):
    """
    Specifies the mapping to curvilinear coordinates.
    
    Takes as input: array_list made by x_coordinates, y_ccordinates in the map 
                    space.

    Returns as output: array_list made by x_coordinates, y_ccordinates in the 
                       physical space.

    Inputs: mC = list composed by two arrays
                 [array ([xc1, xc2, ...]), array([yc1, yc2, ...])]

    Output: pC = list composed by three arrays
                 [array ([xp1, xp2, ...]), array([yp1, yp2, ...]), array([zp1, zp2, ...])]
    """  

    # Radius of the sphere
    r1 = Rsphere
    
    # Number of cell in x and y directions. (x,y) c
    mx = grid.ng[0]
    my = grid.ng[1]

    # Define new list, pC = physical coordinates
    pC = [np.zeros((mx,my))]*3

    for i in range(mx):
        for j in range(my):
            xc = mC[0][i][j]
            yc = mC[1][i][j]

            # Ghost cell values outside of [-3,1]x[-1,1] get mapped to other
            # hemisphere:
            if (xc >= 1.0):
                xc = xc - 4.0
            if (xc <= -3.0):
                xc = xc + 4.0

            if (yc >= 1.0):
                yc = 2.0 - yc
                xc = -2.0 - xc

            if (yc <=-1.0):
                yc = -2.0 - yc
                xc = -2.0 - xc

            if (xc <= -1.0):
                # Points in [-3,-1] map to lower hemisphere - reflect about x=-1
                # to compute x,y mapping and set sgnz appropriately:
                xc = -2.0 - xc
                sgnz = -1.0
            else:
                sgnz = 1.0

            import math
            sgnxc = math.copysign(1.0,xc)
            sgnyc = math.copysign(1.0,yc)

            xc1 = np.abs(xc)
            yc1 = np.abs(yc)
            d = np.maximum(np.maximum(xc1,yc1), 1.0e-10)     

            DD = r1*d*(2.0 - d) / np.sqrt(2.0)
            R = r1
            center = DD - np.sqrt(np.maximum(R**2 - DD**2, 0.0))
            
            xp = DD/d * xc1
            yp = DD/d * yc1

            if (yc1 >= xc1):
                yp = center + np.sqrt(np.maximum(R**2 - xp**2, 0.0))
            else:
                xp = center + np.sqrt(np.maximum(R**2 - yp**2, 0.0))

            # Compute physical coordinates
            pC[2][i][j] = np.sqrt(np.maximum(r1**2 - (xp**2 + yp**2), 0.0))
            pC[0][i][j] = xp*sgnxc
            pC[1][i][j] = yp*sgnyc
            pC[2][i][j] = pC[2][i][j]*sgnz
                  

def qinit(state):
    r"""
    Initialize data with with a Gaussian pulse centered at (x0,y0) with radius 
    r0.
    """
    TO BE DEFINED. 

def shallow_sphere(use_petsc=False,iplot=0,htmlplot=False,outdir='./_output',solver_type='classic'):
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
    
    #===========================================================================
    # Initialize grid and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    # Grid:
    xlower = 0.0
    xupper = 1.0
    mx = 10

    ylower = 0.0
    yupper = 1.0
    my = 5

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])
    dx = grid.d[0]
    dy = grid.d[1]

    # Override default mapc2p function
    grid.mapc2p = mapc2p_sphere
    
    # Define state object
    meqn = 3  # Number of equations
    maux = 16 # Number of auxiliary variables
    state = pyclaw.State(grid,meqn,maux)

    # Set auxiliary variables
    state.aux = init.setaux(mx,my,mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

    # Set initial condition for q
    state.q = init.qinit(mx,my,meqn,mbc,mx,my,xlower,ylower,dx,dy,state.q,maux,aux)
     

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow_sphere)
    print 'Error: ',output







