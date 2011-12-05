#!/usr/bin/env python
# encoding: utf-8

#===========================================================================
# Import libraries
#===========================================================================
import numpy as np


def advection_implicitSharpClawBE(use_petsc=True,iplot=False,htmlplot=False,solver_type='sharpclawBE',outdir='./_output'):
    """
    Example python script for solving the 1d advection equation using the 
    fifth-order WENO scheme and the backward Euler method for time stepping.
    """
    #===========================================================================
    # Import libraries
    #===========================================================================
    # Import only petclaw because the implicit LW always uses SNES to solve the
    # nonlinear algebraic system
    if use_petsc:
        import petclaw as pyclaw
    else:
        raise Exception('Implicit SharpClaw solvers require PETSc!!')

    #===========================================================================
    # Setup solver and solver parameters
    #=========================================================================== 
    solver = pyclaw.ImplicitSharpClawSolverfsolve1D()

    solver.kernel_language = 'Fortran'
    
    from riemann import rp_advection
    solver.mwaves = rp_advection.mwaves
 
    solver.mthbc_lower[0] = pyclaw.BC.periodic
    solver.mthbc_upper[0] = pyclaw.BC.periodic

    solver.time_integrator ='BEuler'

    solver.cfl_desired=1
    solver.cfl_max=1.1
    
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
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    import sys
    import petsc4py
    petclaw_args = [arg for arg in sys.argv[1:] if '=' in arg]
    petsc_args = [arg for arg in sys.argv[1:] if '=' not in arg]
    petsc4py.init(petsc_args)
    petclaw_args.insert(0,sys.argv[0])
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(petclaw_args)
        advection_implicitSharpClawBE(*args,**kwargs)
    else: advection_implicitSharpClawBE()
