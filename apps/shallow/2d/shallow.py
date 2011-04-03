#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid):
    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)

    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    q[0,:,:] = 1.0
    q[1,:,:] = 0.0
    q[2,:,:] = 0.0
    grid.q=q

def shallow2D(iplot=False,petscPlot=False,useCcontroller=True,htmlplot=False)
    """
    Example python script for solving the 2d shallow water equations.
    """

    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from petclaw.solution import Solution
    from petclaw.evolve.clawpack import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

    # Create grid
    mx=100; my=100
    x = Dimension('x',-1.0,1.0,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',-1.0,1.0,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])

    # Define number of equations and BC
    grid.meqn = 3
    grid.mbc = 2

    # Define parameters for simulation
    #g = 1.0 # gravity
    #grid.aux_global['g'] = g
    #from dimsp2 import cparam
    #for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    # Set initial condition (initial solution)
    qinit(grid)
    initial_solution = Solution(grid)


    # Define solver and solver's parameters
    solver = PetClawSolver2D()
    solver.order = 2
    solver.order_trans = 2
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 3
    solver.mathlim = [3]

    # Define controller and controller's parameters
    claw = Controller()
    claw.keep_copy = True
    claw.output_format = 'petsc' # The output format MUST be set to petsc!
    tfinal = 0.1
    claw.tfinal = tfinal
    claw.solution['n'] = initial_solution
    claw.solver = solver

    # Solve problem
    status = claw.run()

    if htmlplot:  plot.plotHTML()
    if petscPlot: plot.plotPetsc(claw)
    if iplot:     plot.plotInteractive()


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        shallow2D(*args,**kwargs)
    else: shallow2D()



    






