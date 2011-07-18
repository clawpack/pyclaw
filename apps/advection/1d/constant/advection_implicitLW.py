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
 
    solver.mthbc_lower[0] = 2
    solver.mthbc_upper[0] = 2

    
    #===========================================================================
    # Initialize grid and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,100)
    grid = pyclaw.Grid(x)
    state = pyclaw.State(grid)
    state.aux_global['u']=1.
    state.meqn = rp_advection.meqn

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
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        error=advection_implicitLW(*args,**kwargs)
        print 'Error: ',error
    else: advection_implicitLW()
