#!/usr/bin/env python
# encoding: utf-8

"""
1D shallow water equations.
"""

    
def shallow1D(iplot=True,petscPlot=False,useController=True,htmlplot=False):
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np
    from petsc4py import PETSc   
    from petclaw.grid import Grid
    from petclaw.grid import Dimension
    from pyclaw.solution import Solution
    from petclaw.evolve.clawpack import PetClawSolver1D
    from pyclaw.controller import Controller
    from petclaw import plot

    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================

    # Grid:
    xlower = -5.0
    xupper = 5.0
    mx = 500
    x = Dimension('x',xlower,xupper,mx,mthbc_lower=3,mthbc_upper=1)
    grid = Grid(x)

    grid.meqn = 2
    grid.t = 0.0

    # Parameters
    grid.aux_global['grav'] = 1.0
    from step1 import cparam
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)


    # Initial solution:
    grid.init_q_petsc_structures() # This must be called before grid.x.center and such can be accessed.
    xCenter = grid.x.center
    q = np.zeros([grid.meqn, len(xCenter)], order = 'F')

    radDam = 0.0
    hl = 3.
    ul = 0.
    hr = 1.
    ur = 0.

    q[0,:] = hl * (grid.p_center[0] <= radDam) + hr * (grid.p_center[0] > radDam)
    q[1,:] = hl*ul * (grid.p_center[0] <= radDam) + hr*ur * (grid.p_center[0] > radDam)
    grid.q=q

    init_solution = Solution(grid)

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    kernelsType = 'F'
    solver = PetClawSolver1D(kernelsType = kernelsType)
    solver.mwaves = 2
    if kernelsType =='P': solver.set_riemann_solver('shallow_roe')
    solver.mthlim = [4]*solver.mwaves

    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    claw = Controller()
    claw.keep_copy = True
    claw.output_format = 'petsc' # The output format MUST be set to petsc!!
    claw.tfinal = 5.0
    claw.solutions['n'] = init_solution
    claw.solver = solver


    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    if htmlplot: plot.plotHTML()
    if petscPlot: plot.plotPetsc(output_object)
    if iplot: plot.plotInteractive()


if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=shallow1D(*args,**kwargs)
    print 'Error: ',error
    

   

