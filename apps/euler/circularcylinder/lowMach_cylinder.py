#!/usr/bin/env python
# encoding: utf-8

#===========================================================================
# Import libraries
#===========================================================================
import numpy as np
from petsc4py import PETSc
import petclaw

def mapc2p_cylinder(grid,mC):
    """
    Maps a rectangular (computational) domain into a circular (physical) domain.
    This is a typical approach to mesh a cylinder and define the physical domain.
    This approach avoids degenerate elements which appear frequently when a
    cylinder is meshed with a structured grid inside a rectangular region.

    The flow is steady (no viscous terms) and, except for the transient, no wake
    passes through the outflow boundary. Therefore, the shape of the domain 
    in combination with the free stream BC should not pose any problem.

    Polar coordinates are used to map the rectangle into a circular region with 
    a circular hole that represents the cylinder.


    Input: array_list made by x_coordinates, y_ccordinates in the map space
      
           mC = list composed by two array [array ([xc1, xc2, ...]), array([yc1, yc2, ...])]



    Output: array_list made by x_coordinates, y_ccordinates in the physical space

            pC = list composed by two array [array ([xp1, xp2, ...]), array([yp1, yp2, ...])]
    """  

    # Polar coordinates (x coordinate = radius, y coordinate = theta)
    nbrCells = len(mC[0])

    # Define new empty list
    pC = []

    # Populate it with the physical coordinates 
    pC.append(mC[0][:]*np.cos(mC[1][:]))
    pC.append(mC[0][:]*np.sin(mC[1][:]))
    
    return pC

def qinit(state,gamma,gamma1):
    """
    Set initial solution.

    Uniform horizontal flow equals to the free stream condition is imposed. 
    Therefore, the initial velocity field close to the cylinder has a normal 
    component to the wall. This will be "corrected" after few pseudo-time steps.
    """
    # Free stream values
    ####################
    # Density 
    rhoinf = 1.0

    # Speed of sound
    cinf = 1.0

    # Mach number
    Minf = 0.1
    
    # Air constant
    R = 287.06

    # Set conserved variables which represents the solution variables of
    # our problem.
    ########################################################################
    # Temperature (no conserved variable)
    Tinf = cinf**2/(gamma*R)

    # Pressure (no conserved variable)
    pinf = rhoinf*cinf**2/gamma

    # Velocity (no conserved variable)
    uinf = Minf*cinf

    # Momentum
    minf = rhoinf*uinf

    # Specific energy
    einf = pinf/(gamma1) + 1/2*rhoinf*uinf**2

    state.q[0,:,:] = rhoinf
    state.q[1,:,:] = mnf
    state.q[2,:,:] = 0.   # No velocity component in the vertical direction!
    state.q[3,:,:] = einf



def lowMach_cylinder(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='classic')
    """
    Run Euler flow over a circular cylinder.
    """

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

    solver.mthbc_lower[0] = pyclaw.BC.reflecting
    solver.mthbc_upper[0] = pyclaw.BC.custom      # Free stream BC (outer circle is very far from the cylinder)
    solver.user_bc_lower  = freestream_lower

    # These 2 boundary segments with the current mapping become the two ends of
    # a C-mesh. It would be nice to have the options two impose two different 
    # kind of BC to the two segments. This is not possible because user_bc_lower
    # has been already taken above. This need to be discussed with the developers.
    solver.mthbc_lower[1] = pyclaw.BC.reflecting  
    solver.mthbc_upper[1] = pyclaw.BC.reflecting


    # The imposition of the BC on the aux array depends on how we want to
    # implement the construction of aux for the Euler equations on map grids!
    solver.mthauxbc_lower[0] =
    solver.mthauxbc_upper[0] = 
    solver.mthauxbc_lower[1] = 
    solver.mthauxbc_upper[1] = 


    solver.mwaves = 3

    solver.dim_split = 1
    solver.order_trans = 2
    solver.order = 2

    solver.dt_initial = 0.1
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.2


    #===========================================================================
    # Initialize grid and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    # Grid:
    xlower = 0.0
    xupper = 1.0
    mx = 40

    ylower = 0.0
    yupper = np.pi*2.0
    my = 120

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])
    grid.mapc2p = mapc2p_annulus # Override default_mapc2p function implemented in grid.py

    meqn = 4  # Number of equations
    state = pyclaw.State(grid,meqn)

    # Set auxiliary gloabl values that will be used in the Riemann solver
    # Ratio cv/cr
    gamma = 1.4
    gamma1 = gamma - 1

    state.aux_global['gamma']= gamma
    state.aux_global['gamma1']= gamma1

    
    # Set initial solution
    # ====================
    qinit(state,gamma,gamma1) # This function is defined above

    
    #
    # SETAUX!!!!
    #


    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.outstyle = 1
    claw.nout = 25
    claw.tfinal = 1.0
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
        error=lowMach_cylinder(*args,**kwargs)
        print 'Error: ',error
    else: lowMach_cylinder()




   


