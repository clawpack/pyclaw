#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(use_petsc=False,kernel_language='Fortran',solver_type='classic',iplot=False,htmlplot=False,outdir='./_output'):
    """
    This example solves the 1-dimensional acoustics equations in a homogeneous
    medium.
    """
    from numpy import sqrt, exp, cos

    #=================================================================
    # Import the appropriate classes, depending on the options passed
    #=================================================================
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else: raise Exception('Unrecognized value of solver_type.')

    #========================================================================
    # Instantiate the solver and define the system of equations to be solved
    #========================================================================
    solver.kernel_language=kernel_language
    from riemann import rp_acoustics
    solver.mwaves=rp_acoustics.mwaves
    if kernel_language=='Python': 
        solver.rp = rp_acoustics.rp_acoustics_1d
 
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.mthbc_lower[0] = pyclaw.BC.periodic
    solver.mthbc_upper[0] = pyclaw.BC.periodic

    #========================================================================
    # Instantiate the grid and set the boundary conditions
    #========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,100)
    grid = pyclaw.Grid(x)
    state = pyclaw.State(grid)

    #========================================================================
    # Set problem-specific variables
    #========================================================================
    rho = 1.0
    bulk = 1.0
    state.aux_global['rho']=rho
    state.aux_global['bulk']=bulk
    state.aux_global['zz']=sqrt(rho*bulk)
    state.aux_global['cc']=sqrt(rho/bulk)

    #========================================================================
    # Set the initial condition
    #========================================================================
    state.meqn=rp_acoustics.meqn
    xc=grid.x.center
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = exp(-beta * (xc-x0)**2) * cos(gamma * (xc - x0))
    state.q[1,:] = 0.
    
    #========================================================================
    # Set up the controller object
    #========================================================================
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = 1.0

    # Solve
    status = claw.run()

    # Plot results
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)
