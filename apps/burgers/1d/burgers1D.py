#!/usr/bin/env python
# encoding: utf-8

def burgers(kernel_language='Fortran',iplot=False,petscPlot=False,useController=True,htmlplot=True):
    """
    Example python script for solving the 1d Burgers equation.
    """

    import numpy as np
    from petsc4py import PETSc
    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.clawpack import PetClawSolver1D
    from pyclaw.controller import Controller
    import math
    from petclaw import plot

    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================
    
    # Grid: 
    x = Dimension('x',0.0,1.0,500,mthbc_lower=2,mthbc_upper=2)
    grid = Grid(x)


    grid.meqn = 1
    grid.t = 0.0

    # Initial solution
    grid.init_q_petsc_structures()
    xc=grid.x.center
    q=np.zeros([len(xc),grid.meqn], order = 'F')
    q[:,0] = np.sin(np.pi*2*xc) + 0.50
    grid.q=q
    grid.aux_global['efix']=True

    initial_solution = Solution(grid)

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    solver = PetClawSolver1D()
    solver.kernel_language = kernel_language
    solver.mwaves = 1
    solver.mthlim = [3]
    if kernel_language=='Python': solver.set_riemann_solver('burgers')



    #===========================================================================
    # Setup controller and controller paramters. Then solve the problem
    #===========================================================================
    if useController:
        claw = Controller()
        claw.keep_copy = True
        # The output format MUST be set to petsc!
        claw.output_format = 'petsc'
        claw.tfinal =0.5
        claw.solutions['n'] = initial_solution
        claw.solver = solver

        status = claw.run()
        output_object=claw

    else:
        sol = {"n":init_solution}
        solver.evolve_to_time(sol,0.25)
        output_object=sol["n"]

    if petscPlot:
        plot.plotPetsc(output_object)

    if iplot:
        plot.plotInteractive()

    if htmlplot:
        plot.plotHTML()


    return output_object



if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=burgers(*args,**kwargs)
    print 'Error: ',error

