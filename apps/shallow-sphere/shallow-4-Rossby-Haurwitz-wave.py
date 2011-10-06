#!/usr/bin/env python
# encoding: utf-8

"""
2D shallow water equations on a spherical surface.
The approximation of the three-dimensional equations
is restricted to the surface of the sphere.
"""

import numpy as np
from petclaw import plot
import matplotlib.pyplot as plt

# Nondimensionalized radius of the earth
Rsphere = 1.0

def fortran_src_wrapper(solver,solution,t,dt):
    """
    Wraps Fortran src2.f routine.
    """

    grid = solution.grid
    state = solution.state

    mx,my = grid.ng[0], grid.ng[1]
    meqn = state.meqn
    mbc = solver.mbc
    xlower,ylower = grid.lower[0], grid.lower[1]
    dx,dy = grid.d[0],grid.d[1]
    q = state.q
    maux = state.maux
    aux = state.aux

    import problem
    state.q=problem.src2(mx,my,mbc,xlower,ylower,dx,dy,q,aux,t,dt,Rsphere)


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
    """  

    # Radius of the sphere
    r1 = Rsphere
    
    # Number of cell in x and y directions. (x,y) c
    mx = grid.ng[0]
    my = grid.ng[1]

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
            zp = np.sqrt(np.maximum(r1**2 - (xp**2 + yp**2), 0.0))
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
    """

    # Nondimensionalized radius of the earth
    r1 = Rsphere

    # Number of cell in x and y directions. (x,y) c
    mx = grid.ng[0]
    my = grid.ng[1]

    sgnz = np.ones((mx,my))

    # 2D coordinates in the computational domain
    ############################################
    xc = mC[0][:][:]
    yc = mC[1][:][:]


    # 3D coordinates in the physical domain
    #######################################
    # yc < -1 => second copy of sphere:
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
    D = r1*d*(2-d) / np.sqrt(2)
    R = r1*np.ones((np.shape(d)))

    center = D - np.sqrt(R**2 - D**2)
    xp = D/d * xc1
    yp = D/d * yc1

    ij = np.where(yc1==d)
    yp[ij] = center[ij] + np.sqrt(R[ij]**2 - xp[ij]**2)
    ij = np.where(xc1==d)
    xp[ij] = center[ij] + np.sqrt(R[ij]**2 - yp[ij]**2)
    
    # Define new list of numpy array, pC = physical coordinates
    pC = []

    xp = np.sign(xc) * xp
    yp = np.sign(yc) * yp
    zp = sgnz * np.sqrt(r1**2 - (xp**2 + yp**2))
    
    pC.append(xp)
    pC.append(yp)
    pC.append(zp)

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
    h0 = 8.e3        
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


def qbc_lower_y(grid,dim,t,qbc,mbc):
    """
    Impose periodic boundary condition to q at the bottom boundary for the 
    sphere. This function is vectorized.
    """
    if dim == grid.dimensions[1]:
        for j in range(mbc):
            qbc1D = np.copy(qbc[:,:,2*mbc-1-j])
            qbc[:,:,j] = qbc1D[:,::-1]
    else:
        raise Exception('Aux custum BC for this boundary is not appropriate!')



def qbc_upper_y(grid,dim,t,qbc,mbc):
    """
    Impose periodic boundary condition to q at the top boundary for the sphere.
    This function is vectorized.
    """
    if dim == grid.dimensions[1]:
        my = grid.ng[1]
        for j in range(mbc):
            qbc1D = np.copy(qbc[:,:,my+mbc-1-j])
            qbc[:,:,my+mbc+j] = qbc1D[:,::-1]
    else:
        raise Exception('Custum BC for this boundary is not appropriate!')


def auxbc_lower_y(grid,dim,t,auxbc,mbc):
    """
    Impose periodic boundary condition to aux at the bottom boundary for the 
    sphere. This function is vectorized.
    """
    import problem
    
    mx,my = grid.ng[0], grid.ng[1]
    xlower,ylower = grid.lower[0], grid.lower[1]
    dx,dy = grid.d[0],grid.d[1]

    auxbc = problem.setaux(mx,my,mbc,mx,my,xlower,ylower,dx,dy,auxbc,Rsphere)

    #if dim == grid.dimensions[1]:
    #    for j in range(mbc):
    #        auxbc1D = auxbc[:,:,2*mbc-1-j].copy()
    #        auxbc[:,:,j] = auxbc1D[:,::-1]
    #else:
    #    raise Exception('Aux custum BC for this boundary is not appropriate!')


def auxbc_upper_y(grid,dim,t,auxbc,mbc):
    """
    Impose periodic boundary condition to aux at the top boundary for the 
    sphere. This function is vectorized.
    """
    import problem

    mx,my = grid.ng[0], grid.ng[1]
    xlower,ylower = grid.lower[0], grid.lower[1]
    dx,dy = grid.d[0],grid.d[1]

    auxbc = problem.setaux(mx,my,mbc,mx,my,xlower,ylower,dx,dy,auxbc,Rsphere)

    #if dim == grid.dimensions[1]:
    #    my = grid.ng[1]
    #    for j in range(mbc):
    #        auxbc1D = auxbc[:,:,my+mbc-1-j].copy()
    #        auxbc[:,:,my+mbc+j] = auxbc1D[:,::-1]
    #else:
    #    raise Exception('Aux custum BC for this boundary is not appropriate!')


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

    # Set boundary conditions
    # =======================
    solver.mthbc_lower[0] = pyclaw.BC.periodic
    solver.mthbc_upper[0] = pyclaw.BC.periodic
    solver.mthbc_lower[1] = pyclaw.BC.custom
    solver.mthbc_upper[1] = pyclaw.BC.custom

    solver.user_bc_lower = qbc_lower_y
    solver.user_bc_upper = qbc_upper_y

    solver.mthauxbc_lower[0] = pyclaw.BC.periodic
    solver.mthauxbc_upper[0] = pyclaw.BC.periodic
    solver.mthauxbc_lower[1] = pyclaw.BC.custom
    solver.mthauxbc_upper[1] = pyclaw.BC.custom

    solver.user_aux_bc_lower = auxbc_lower_y
    solver.user_aux_bc_upper = auxbc_upper_y


    # Order of the solver.
    # ====================
    solver.order = 2

    # Dimensional splitting ?
    # =======================
    solver.dim_split = 0
 
    # Transverse increment waves and transverse correction waves are computed 
    # and propagated.
    # =======================================================================
    solver.order_trans = 2
    
    # Number of waves in each Riemann solution
    # ========================================
    solver.mwaves = 3

    # Use source splitting method
    # ===========================
    solver.src_split = 2

    # Set variable dt
    # ===============
    solver.dt_variable = 1


    # Set source function
    # ===================
    solver.src = fortran_src_wrapper

    # Set time stepping parameters
    # ============================
    solver.dt_initial = 0.1
    solver.cfl_max = 1.0
    solver.cfl_desired = 0.9

    # Set the same limiter for the mwaves waves. 
    # ==========================================
    solver.limiters = pyclaw.limiters.tvd.MC


    #===========================================================================
    # Initialize grid and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    # Grid:
    # ====
    xlower = -3.0
    xupper = 1.0
    mx = 200

    ylower = -1.0
    yupper = 1.0
    my = 100

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])
    dx = grid.d[0]
    dy = grid.d[1]

    import classic2
    classic2.comxyt.dxcom = dx
    classic2.comxyt.dycom = dy
    classic2.sw.g = 11489.57219

    # Override default mapc2p function
    # ================================
    grid.mapc2p = mapc2p_sphere_vectorized
        
    # Define state object
    # ===================
    meqn = 4  # Number of equations
    maux = 16 # Number of auxiliary variables
    state = pyclaw.State(grid,meqn,maux)

    # Set auxiliary variables
    # =======================
    import problem

    # 1) Call to simplified Fortran function
    #state.aux = problem.setaux(xlower,ylower,dx,dy,state.aux,Rsphere,mx,my,maux)

    # 2) Call to original Fortran function
    # TO USE THIS ONE: RECNAME qinitOrig.f to qinit.f and recompile (make)
    # THIS OPTION WILL BE REMOVED SOON.
    mbc = 2
    auxtmp = np.ndarray(shape=(maux,mx+2*mbc,my+2*mbc), dtype=float, order='F')
    auxtmp = problem.setaux(mx,my,mbc,mx,my,xlower,ylower,dx,dy,auxtmp,Rsphere)
    state.aux[:,:,:] = auxtmp[:,mbc:-mbc,mbc:-mbc]

    # Set index for capa
    state.mcapa = 0

    # Set initial condition for q
    #############################
    # 1) Call to simplified Fortran function
    #state.q = problem.qinit(xlower,ylower,dx,dy,state.q,state.aux,Rsphere,meqn,mx,my,maux)

    
    # 2) Call to original Fortran function
    # TO USE THIS ONE: RECNAME qinitOrig.f to qinit.f and recompile (make)
    # THIS OPTION WILL BE REMOVED SOON.
    qtmp = np.ndarray(shape=(meqn,mx+2*mbc,my+2*mbc), dtype=float, order='F')
    qtmp = problem.qinit(mx,my,mbc,mx,my,xlower,ylower,dx,dy,qtmp,auxtmp,Rsphere)
    state.q[:,:,:] = qtmp[:,mbc:-mbc,mbc:-mbc]

    # 3) call to python function define above
    #qinit(state,mx,my)


    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.outstyle = 1
    claw.nout = 24
    claw.tfinal = 24
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
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow_sphere)
    print 'Error: ',output







