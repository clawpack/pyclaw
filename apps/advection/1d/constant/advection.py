#!/usr/bin/env python
# encoding: utf-8
def advection(kernelsType='P',iplot=True,petscPlot=False,useController=True):
    """
    Example python script for solving the 1d advection equation.
    """

    import numpy as np
    from petsc4py import PETSc
    from petclaw.grid import Dimension, Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.solver import PetClawSolver1D
    from pyclaw.controller import Controller
    from petclaw import plot

    # Initialize grids and solutions
    x = Dimension('x',0.0,1.0,100,mthbc_lower=2,mthbc_upper=2)
    grid = Grid(x)
    grid.aux_global['u']=-1.
    if kernelsType=='F': 
        from step1 import comrp
        comrp.u = grid.aux_global['u']
    grid.meqn = 1
    grid.t = 0.0
    grid.init_q_petsc_structures()

    xc=grid.x.center
    beta=100; gamma=0; x0=0.75
    grid.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

    initial_solution = Solution(grid)

    # Solver setup
    solver = PetClawSolver1D(kernelsType = 'P')
    solver.mwaves = 1
    if kernelsType=='P': solver.set_riemann_solver('advection')
    solver.mthlim=[4]

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

    return output_object

if __name__=="__main__":
    advection()
