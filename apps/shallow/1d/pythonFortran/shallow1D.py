#!/usr/bin/env python
# encoding: utf-8

"""
1D shallow water equations.
"""

    
def shallow1D(use_PETSc=False,kernel_language='Fortran',iplot=True,userController=True,petscPlot=False,htmlplot=False,outdir='./_output'):
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
    x = Dimension('x',xlower,xupper,mx,mthbc_lower=1,mthbc_upper=1)
    grid = Grid(x)

    grid.meqn = 2
    grid.t = 0.0

    # Parameters
    grid.aux_global['grav'] = 1.0
    from classic1 import cparam
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)


    # Initial solution:
    grid.init_q_petsc_structures() # This must be called before grid.x.center and such can be accessed.
    xCenter = grid.x.center
    q = np.zeros([grid.meqn, len(xCenter)], order = 'F')

    damRadius = 0.0
    hl = 3.
    ul = 0.
    hr = 1.
    ur = 0.

    q[0,:] = hl * (grid.p_center[0] <= damRadius) + hr * (grid.p_center[0] > damRadius)
    q[1,:] = hl*ul * (grid.p_center[0] <= damRadius) + hr*ur * (grid.p_center[0] > damRadius)
    grid.q=q

    init_solution = Solution(grid)

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    solver = PetClawSolver1D()
    solver.mwaves = 2
    solver.mthlim = [4]*solver.mwaves
    solver.kernel_language=kernel_language
    if kernel_language =='Python': 
        solver.set_riemann_solver('shallow_roe')
        grid.aux_global['g'] = 1.0
        grid.aux_global['efix'] = False


    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    claw = Controller()
    claw.keep_copy = True
    claw.tfinal = 2.0
    claw.solutions['n'] = init_solution
    claw.solver = solver
    claw.outdir = outdir


    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    #===========================================================================
    # Plot results
    #===========================================================================
    if iplot:     plot.plotInteractive(outdir=outdir,format=claw.output_format)
    if htmlplot:  plot.plotHTML(outdir=outdir,format=claw.output_format)


if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=shallow1D(*args,**kwargs)
    print 'Error: ',error
    

   

