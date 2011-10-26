#!/usr/bin/env python
# encoding: utf-8

#===========================================================================
# Import libraries
#===========================================================================
import numpy as np


def advection_implicitLW(use_petsc=True,iplot=False,htmlplot=False,solver_type='classic',outdir='./_output'):
    """
    Example python script for solving the 1d advection equation using the 
    implicit Lax-Wendroff scheme.
    """
    #===========================================================================
    # Import libraries
    #===========================================================================
    # Import only petclaw because the implicit LW always uses SNES to solve the
    # nonlinear algebraic system
    if use_petsc:
        import petclaw as pyclaw
    else:
        raise Exception('The Implicit LW method requires PETSc!!')

    #===========================================================================
    # Setup solver and solver parameters
    #=========================================================================== 
    solver = pyclaw.ImplicitClawSolver1D()

    solver.kernel_language = 'Fortran'
    
    from riemann import rp_advection
    solver.mwaves = rp_advection.mwaves
 
    solver.bc_lower[0] = 2
    solver.bc_upper[0] = 2

    
    #===========================================================================
    # Initialize grid and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,100)
    grid = pyclaw.Grid(x)
    meqn = 1
    state = pyclaw.State(grid,meqn)
    state.aux_global['u']=1.

    xc = grid.x.center
    beta = 100
    gamma = 0
    x0 = 0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir

    claw.tfinal = 1.0

    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    #===========================================================================
    # Plot results
    #===========================================================================
    #if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    #if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":

    from pyclaw.util import run_app_from_main
    output = run_app_from_main(advection_implicitLW)
    print 'Error: ',output

