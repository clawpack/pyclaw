#!/usr/bin/env python
# encoding: utf-8

"""
2D shallow water equations on a spherical surface. The approximation of the 
three-dimensional equations is restricted to the surface of the sphere. 
Therefore only the solution on the surface is updated. 

Reference: Logically Rectangular Grids and Finite Volume Methods for PDEs in 
           Circular and Spherical Domains. 
           By Donna A. Calhoun, Christiane Helzel, and Randall J. LeVeque
           SIAM Review 50 (2008), 723-752. 
"""

# Import library
# ==============
import numpy as np

from numpy import *
import os

try:
    import problem
    import classic2
except ImportError:
    import warnings
    warnings.warn("missing extension modules, running python setup.py build_ext -i")
    import subprocess
    subprocess.check_call('python setup.py build_ext -i', shell=True)
    warnings.warn("missing extension modules built by running python setup.py build_ext -i")
    try:
        import problem
        import classic2
    except ImportError:
        import sys
        print >> sys.stderr, "***\nUnable to import problem module or automatically build, try running (in the directory of this file):\n python setup.py build_ext -i\n***"
        raise


# Parameters used by the following routines
# =========================================

# Nondimensionalized radius of the earth
Rsphere = 1.0


def fortran_src_wrapper(solver,state,dt):
    """
    Wraps Fortran src2.f routine. 
    src2.f contains the discretization of the source term.
    """
    # Some simplifications
    grid = state.grid

    # Get parameters and variables that have to be passed to the fortran src2
    # routine.
    mx, my = grid.num_cells[0], grid.num_cells[1]
    num_ghost = solver.num_ghost
    xlower, ylower = grid.lower[0], grid.lower[1]
    dx, dy = grid.delta[0], grid.delta[1]
    q = state.q
    aux = state.aux
    t = state.t

    # Call src2 function
    state.q = problem.src2(mx,my,num_ghost,xlower,ylower,dx,dy,q,aux,t,dt,Rsphere)


def mapc2p_sphere_nonvectorized(grid,mC):
    """
    Maps to points on a sphere of radius Rsphere. Nonvectorized version (slow).
    
    Takes as input: array_list made by x_coordinates, y_ccordinates in the map 
                    space.

    Returns as output: array_list made by x_coordinates, y_ccordinates in the 
                       physical space.

    Inputs: mC = list composed by two arrays
                 [array ([xc1, xc2, ...]), array([yc1, yc2, ...])]

    Output: pC = list composed by three arrays
                 [array ([xp1, xp2, ...]), array([yp1, yp2, ...]), array([zp1, zp2, ...])]

    NOTE: this function is not used in the standard script.
    """  
    # Import library            
    import math

    # Get number of cells in both directions
    mx, my = grid.num_cells[0], grid.num_cells[1]

    # Define new list of numpy array, pC = physical coordinates
    pC = []

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

            sgnxc = math.copysign(1.0,xc)
            sgnyc = math.copysign(1.0,yc)

            xc1 = np.abs(xc)
            yc1 = np.abs(yc)
            d = np.maximum(np.maximum(xc1,yc1), 1.0e-10)     

            DD = Rsphere*d*(2.0 - d) / np.sqrt(2.0)
            R = Rsphere
            centers = DD - np.sqrt(np.maximum(R**2 - DD**2, 0.0))
            
            xp = DD/d * xc1
            yp = DD/d * yc1

            if (yc1 >= xc1):
                yp = centers + np.sqrt(np.maximum(R**2 - xp**2, 0.0))
            else:
                xp = centers + np.sqrt(np.maximum(R**2 - yp**2, 0.0))

            # Compute physical coordinates
            zp = np.sqrt(np.maximum(Rsphere**2 - (xp**2 + yp**2), 0.0))
            pC.append(xp*sgnxc)
            pC.append(yp*sgnyc)
            pC.append(zp*sgnz)

    return pC


def mapc2p_sphere_vectorized(grid,mC):
    """
    Maps to points on a sphere of radius Rsphere. Vectorized version (fast).  

    Takes as input: array_list made by x_coordinates, y_ccordinates in the map 
                    space.

    Returns as output: array_list made by x_coordinates, y_ccordinates in the 
                       physical space.

    Inputs: mC = list composed by two arrays
                 [array ([xc1, xc2, ...]), array([yc1, yc2, ...])]

    Output: pC = list composed by three arrays
                 [array ([xp1, xp2, ...]), array([yp1, yp2, ...]), array([zp1, zp2, ...])]

    NOTE: this function is used in the standard script.
    """

    # Get number of cells in both directions
    mx, my = grid.num_cells[0], grid.num_cells[1]
    
    # 2D array useful for the vectorization of the function
    sgnz = np.ones((mx,my))

    # 2D coordinates in the computational domain
    xc = mC[0][:][:]
    yc = mC[1][:][:]

    # Compute 3D coordinates in the physical domain
    # =============================================

    # Note: yc < -1 => second copy of sphere:
    ij2 = np.where(yc < -1.0)
    xc[ij2] = -xc[ij2] - 2.0;
    yc[ij2] = -yc[ij2] - 2.0;

    ij = np.where(xc < -1.0)
    xc[ij] = -2.0 - xc[ij]
    sgnz[ij] = -1.0;
    xc1 = np.abs(xc)
    yc1 = np.abs(yc)
    d = np.maximum(xc1,yc1)
    d = np.maximum(d, 1e-10)
    D = Rsphere*d*(2-d) / np.sqrt(2)
    R = Rsphere*np.ones((np.shape(d)))

    centers = D - np.sqrt(R**2 - D**2)
    xp = D/d * xc1
    yp = D/d * yc1

    ij = np.where(yc1==d)
    yp[ij] = centers[ij] + np.sqrt(R[ij]**2 - xp[ij]**2)
    ij = np.where(xc1==d)
    xp[ij] = centers[ij] + np.sqrt(R[ij]**2 - yp[ij]**2)
    
    # Define new list of numpy array, pC = physical coordinates
    pC = []

    xp = np.sign(xc) * xp
    yp = np.sign(yc) * yp
    zp = sgnz * np.sqrt(Rsphere**2 - (xp**2 + yp**2))
    
    pC.append(xp)
    pC.append(yp)
    pC.append(zp)

    return pC


def qinit(state,mx,my):
    r"""
    Initialize solution with 4-Rossby-Haurwitz wave.

    NOTE: this function is not used in the standard script.
    """
    # Parameters
    a = 6.37122e6     # Radius of the earth
    Omega = 7.292e-5  # Rotation rate
    G = 9.80616       # Gravitational acceleration

    K = 7.848e-6   
    t0 = 86400.0     
    h0 = 8.e3         # Minimum fluid height at the poles        
    R = 4.0

    # Compute the the physical coordinates of the cells' centerss
    state.grid.compute_p_centers(recompute=True)
 
    for i in range(mx):
        for j in range(my):
            xp = state.grid.p_centers[0][i][j]
            yp = state.grid.p_centers[1][i][j]
            zp = state.grid.p_centers[2][i][j]

            rad = np.maximum(np.sqrt(xp**2 + yp**2),1.e-6)

            if (xp >= 0.0 and yp >= 0.0):
                theta = np.arcsin(yp/rad) 
            elif (xp <= 0.0 and yp >= 0.0):
                theta = np.pi - np.arcsin(yp/rad)
            elif (xp <= 0.0 and yp <= 0.0):
                 theta = -np.pi + np.arcsin(-yp/rad)
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
            # ==================================================
            Uin = np.zeros(3)

            # Longitude (angular) velocity component
            Uin[0] = (K*np.cos(yp)+K*np.cos(yp)**(R-1.)*( R*np.sin(yp)**2.0 - \
                     np.cos(yp)**2.0)*np.cos(R*xp))*t0

            # Latitude (angular) velocity component
            Uin[1] = (-K*R*np.cos(yp)**(R-1.0)*np.sin(yp)*np.sin(R*xp))*t0

            # Radial velocity component
            Uin[2] = 0.0 # The fluid does not enter in the sphere
            

            # Calculate velocity vetor in cartesian coordinates
            # =================================================
            Uout = np.zeros(3)

            Uout[0] = (-np.sin(xp)*Uin[0]-np.sin(yp)*np.cos(xp)*Uin[1])
            Uout[1] = (np.cos(xp)*Uin[0]-np.sin(yp)*np.sin(xp)*Uin[1])
            Uout[2] = np.cos(yp)*Uin[1]

            # Set the initial condition
            # =========================
            state.q[0,i,j] =  h0/a + (a/G)*( bigA + bigB*np.cos(R*xp) + \
                              bigC*np.cos(2.0*R*xp))
            state.q[1,i,j] = state.q[0,i,j]*Uout[0] 
            state.q[2,i,j] = state.q[0,i,j]*Uout[1] 
            state.q[3,i,j] = state.q[0,i,j]*Uout[2] 


def qbc_lower_y(state,dim,t,qbc,num_ghost):
    """
    Impose periodic boundary condition to q at the bottom boundary for the 
    sphere. This function does not work in parallel.
    """
    for j in range(num_ghost):
        qbc1D = np.copy(qbc[:,:,2*num_ghost-1-j])
        qbc[:,:,j] = qbc1D[:,::-1]


def qbc_upper_y(state,dim,t,qbc,num_ghost):
    """
    Impose periodic boundary condition to q at the top boundary for the sphere.
    This function does not work in parallel.
    """
    my = state.grid.num_cells[1]
    for j in range(num_ghost):
        qbc1D = np.copy(qbc[:,:,my+num_ghost-1-j])
        qbc[:,:,my+num_ghost+j] = qbc1D[:,::-1]


def auxbc_lower_y(state,dim,t,auxbc,num_ghost):
    """
    Impose periodic boundary condition to aux at the bottom boundary for the 
    sphere.
    """
    grid=state.grid

    # Get parameters and variables that have to be passed to the fortran src2
    # routine.
    mx, my = grid.num_cells[0], grid.num_cells[1]
    xlower, ylower = grid.lower[0], grid.lower[1]
    dx, dy = grid.delta[0],grid.delta[1]

    # Impose BC
    auxtemp = auxbc.copy()
    auxtemp = problem.setaux(mx,my,num_ghost,mx,my,xlower,ylower,dx,dy,auxtemp,Rsphere)
    auxbc[:,:,:num_ghost] = auxtemp[:,:,:num_ghost]

def auxbc_upper_y(state,dim,t,auxbc,num_ghost):
    """
    Impose periodic boundary condition to aux at the top boundary for the 
    sphere. 
    """
    grid=state.grid

    # Get parameters and variables that have to be passed to the fortran src2
    # routine.
    mx, my = grid.num_cells[0], grid.num_cells[1]
    xlower, ylower = grid.lower[0], grid.lower[1]
    dx, dy = grid.delta[0],grid.delta[1]
    
    # Impose BC
    auxtemp = auxbc.copy()
    auxtemp = problem.setaux(mx,my,num_ghost,mx,my,xlower,ylower,dx,dy,auxtemp,Rsphere)
    auxbc[:,:,-num_ghost:] = auxtemp[:,:,-num_ghost:]


def shallow_4_Rossby_Haurwitz(use_petsc=False,solver_type='classic',iplot=0,htmlplot=False,outdir='./_output', disable_output=False):

    # Import pyclaw module
    if use_petsc:
        raise Exception("petclaw does not currently support mapped grids (go bug Lisandro who promised to implement them)")

    if solver_type != 'classic':
        raise Exception("Only Classic-style solvers (solver_type='classic') are supported")

    from clawpack import pyclaw

    #===========================================================================
    # Set up solver and solver parameters
    #===========================================================================
    solver = pyclaw.ClawSolver2D()
    from clawpack import riemann
    solver.rp = riemann.rp2_shallow_sphere
    solver.fmod = classic2

    # Set boundary conditions
    # =======================

    # Conserved variables
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.bc_lower[1] = pyclaw.BC.custom  # Custom BC for sphere
    solver.bc_upper[1] = pyclaw.BC.custom  # Custom BC for sphere

    solver.user_bc_lower = qbc_lower_y
    solver.user_bc_upper = qbc_upper_y

    # Auxiliary array
    solver.aux_bc_lower[0] = pyclaw.BC.periodic
    solver.aux_bc_upper[0] = pyclaw.BC.periodic
    solver.aux_bc_lower[1] = pyclaw.BC.custom  # Custom BC for sphere
    solver.aux_bc_upper[1] = pyclaw.BC.custom  # Custom BC for sphere

    solver.user_aux_bc_lower = auxbc_lower_y
    solver.user_aux_bc_upper = auxbc_upper_y


    # Dimensional splitting ?
    # =======================
    solver.dimensional_split = 0
 
    # Transverse increment waves and transverse correction waves are computed 
    # and propagated.
    # =======================================================================
    solver.transverse_waves = 2
    
    # Number of waves in each Riemann solution
    # ========================================
    solver.num_waves = 3

    # Use source splitting method
    # ===========================
    solver.source_split = 2

    # Set source function
    # ===================
    solver.step_source = fortran_src_wrapper

    # Set the limiter for the waves
    # =============================
    solver.limiters = pyclaw.limiters.tvd.MC


    #===========================================================================
    # Initialize domain and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    # Domain:
    # ====
    xlower = -3.0
    xupper = 1.0
    mx = 40

    ylower = -1.0
    yupper = 1.0
    my = 20

    # Check whether or not the even number of cells are used in in both 
    # directions. If odd numbers are used a message is print at screen and the 
    # simulation is interrputed.
    if(mx % 2 != 0 or my % 2 != 0):
        import sys
        message = 'Please, use even numbers of cells in both direction. ' \
                  'Only even numbers allow to impose correctly the boundary ' \
                  'conditions!'
        print message
        sys.exit(0)


    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    domain = pyclaw.Domain([x,y])
    dx = domain.grid.delta[0]
    dy = domain.grid.delta[1]

    # Define some parameters used in Fortran common blocks 
    solver.fmod.comxyt.dxcom = dx
    solver.fmod.comxyt.dycom = dy
    solver.fmod.sw.g = 11489.57219  
    solver.rp.comxyt.dxcom = dx
    solver.rp.comxyt.dycom = dy
    solver.rp.sw.g = 11489.57219  

    # Define state object
    # ===================
    num_eqn = 4  # Number of equations
    num_aux = 16 # Number of auxiliary variables
    state = pyclaw.State(domain,num_eqn,num_aux)

    # Override default mapc2p function
    # ================================
    state.grid.mapc2p = mapc2p_sphere_vectorized
        

    # Set auxiliary variables
    # =======================
    
    # Get lower left corner coordinates 
    xlower,ylower = state.grid.lower[0],state.grid.lower[1]

    num_ghost = 2
    auxtmp = np.ndarray(shape=(num_aux,mx+2*num_ghost,my+2*num_ghost), dtype=float, order='F')
    auxtmp = problem.setaux(mx,my,num_ghost,mx,my,xlower,ylower,dx,dy,auxtmp,Rsphere)
    state.aux[:,:,:] = auxtmp[:,num_ghost:-num_ghost,num_ghost:-num_ghost]

    # Set index for capa
    state.index_capa = 0

    # Set initial conditions
    # ====================== 
    # 1) Call fortran function
    qtmp = np.ndarray(shape=(num_eqn,mx+2*num_ghost,my+2*num_ghost), dtype=float, order='F')
    qtmp = problem.qinit(mx,my,num_ghost,mx,my,xlower,ylower,dx,dy,qtmp,auxtmp,Rsphere)
    state.q[:,:,:] = qtmp[:,num_ghost:-num_ghost,num_ghost:-num_ghost]

    # 2) call python function define above
    #qinit(state,mx,my)


    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.output_style = 1
    claw.num_output_times = 10
    claw.tfinal = 10
    claw.solution = pyclaw.Solution(state,domain)
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

    return claw

        
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow_4_Rossby_Haurwitz)
