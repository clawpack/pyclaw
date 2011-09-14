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

# Parameters
Rsphere = 1.0

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

    # Define new list of numpy array, pC = physical coordinates
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

            if (yc <= -1.0):
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

    return pC


def qinit(state,mx,my):
    r"""
    Initialize data with 4-Rossby-Haurwitz wave.
    """
    # Parameters
    a = 6.37122e6     # Radius of the earth 
    K = 7.848e-6   
    Omega = 7.292e-5  # Rotation rate
    G = 9.80616       # Gravitational acceleration
    t0 = 86400.0     
    h0 = 8.0e3        
    R = 4.0

    # Compute the the physical coordinates of the cells' centers
    state.grid.compute_p_center(recompute=True)
  
    for i in range(mx):
        for j in range(my):
            xp = state.grid.p_center[0][i][j]
            yp = state.grid.p_center[1][i][j]
            zp = state.grid.p_center[2][i][j]

            rad = np.maximum(np.sqrt(xp**2 + yp**2),1.e-6)

            if (xp >= 0.0 and yp >= 0.0):
                theta = np.arcsin(yp/rad) 
            elif (xp <= 0.0 and yp >= 0.0):
                theta = np.pi - np.arcsin(yp/rad)
            elif (xp <= 0.0 and yp <= 0.0):
                 theta = -pi + np.arcsin(-yp/rad)
            elif (xp >= 0.0 and yp <= 0.0):
                theta = -np.arcsin(-yp/rad)

            # Compute phi, at north pole: pi/2 at south pool: -pi/2
            if (zp >= 0.0): 
                phi =  np.arcsin(zp/Rsphere) 
            else:
                phi = -np.arcsin(-zp/Rsphere)  
        
            xp = theta 
            yp = phi 


            bigA = 0.5*K*(2.0*Omega + K)*np.cos(yp)**2.0 + \
                   0.25*K*K*np.cos(yp)**(2.0*R)*((1.0*R+1.0)*np.cos(yp)**2.0 + \
                   (2.0*R*R - 1.0*R - 2.0) - 2.0*R*R*(np.cos(yp))**(-2.0))
            bigB = (2.0*(Omega + K)*K)/((1.0*R + 1.0)*(1.0*R + 2.0)) * \
                   np.cos(yp)**R*( (1.0*R*R + 2.0*R + 2.0) - \
                   (1.0*R + 1.0)**(2)*np.cos(yp)**2 )
            bigC = 0.25*K*K*np.cos(yp)**(2*R)*( (1.0*R + 1.0)* \
                   np.cos(yp)**2 - (1.0*R + 2.0))


            # Calculate local longitude-latitude velocity vector
            ####################################################
            Uin = np.zeros(3)

            # Longitude (angular) velocity component
            Uin[0] = (K*np.cos(yp)+K*np.cos(yp)**(R-1.)*( R*np.sin(yp)**2.0 - \
                     np.cos(yp)**2.0)*np.cos(R*xp))*t0

            # Latitude (angular) velocity component
            Uin[1] = (-K*R*np.cos(yp)**(R-1.0)*np.sin(yp)*np.sin(R*xp))*t0

            # Radial velocity component
            Uin[2] = 0.0 # The fluid does not enter in the sphere
            

            # Calculate velocity vetor in cartesian coordinates
            ###################################################
            Uout = np.zeros(3)

            Uout[0] = (-np.sin(xp)*Uin[0]-np.sin(yp)*np.cos(xp)*Uin[1])
            Uout[1] = (np.cos(xp)*Uin[0]-np.sin(yp)*np.sin(xp)*Uin[1])
            Uout[2] = np.cos(yp)*Uin[1]

            # Set the initial condition             
            state.q[0,i,j] =  h0/a + (a/G)*( bigA + bigB*np.cos(R*xp) + \
                              bigC*np.cos(2.0*R*xp))
            state.q[1,i,j] = state.q[0,i,j]*Uout[0] 
            state.q[2,i,j] = state.q[0,i,j]*Uout[1] 
            state.q[3,i,j] = state.q[0,i,j]*Uout[2] 

    
    

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
    meqn = 4  # Number of equations
    maux = 16 # Number of auxiliary variables
    state = pyclaw.State(grid,meqn,maux)

    # Set auxiliary variables
    #########################
    import init
    #mbc = 2 # This is not very good because the user should not worry about the 
            # number of BC (which are solver dependent) 
    state.aux = init.setaux(mx,my,xlower,ylower,dx,dy,state.aux,Rsphere)
    #state.aux[:,:,:] = 0.0
    #print state.aux

    # Set initial condition for q
    #############################

    # 1) call to Fortran function
    #state.q = init.qinit(mx,my,xlower,ylower,dx,dy,state.q,state.aux,Rsphere)
        
    # 2) call to python function define above
    qinit(state,mx,my)
    print state.q


     

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow_sphere)
    print 'Error: ',output







