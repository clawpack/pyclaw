#!/usr/bin/env python
# encoding: utf-8

#===========================================================================
# Import libraries
#===========================================================================
import numpy as np
from petsc4py import PETSc



def qinit(state):

    # Set initial conditions for q.
    # Sample scalar equation with data that is piecewise constant with
    # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    #     0.1  otherwise
    
    x = state.grid.x.center
    y = state.grid.y.center
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] > 0.0 and x[i] < 0.5 and y[j]>0.0 and y[j] < 0.5:
                state.q[:,i,j] = 1.0
            else:
                state.q[:,i,j] = 0.1
                
def advection2D(iplot=False,use_petsc=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d advection equation.
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
    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()


    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic

    solver.mwaves = 1

    solver.dim_split = 0

    solver.cfl_max=1.0
    solver.cfl_desired = 0.9

    solver.mthlim = pyclaw.limiters.tvd.vanleer


    #===========================================================================
    # Initialize grids, then initialize the solution associated to the grid and
    # finally initialize aux array
    #===========================================================================

    # Grid:
    mx=50; my=50
    x = pyclaw.Dimension('x',0.0,1.0,mx)
    y = pyclaw.Dimension('y',0.0,1.0,my)
    grid = pyclaw.Grid([x,y])

    meqn = 1
    state = pyclaw.State(grid,meqn)

    state.aux_global['u'] = 0.5 # Parameters (global auxiliary variables)
    state.aux_global['v'] = 1.0

    # Initial solution
    # ================
    qinit(state) # This function is defined above


    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 2.0
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
    output = run_app_from_main(advection2D)
