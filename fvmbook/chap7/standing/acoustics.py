#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(use_petsc=False,kernel_language='Fortran',solver_type='classic',iplot=False,htmlplot=False,outdir='./_output',weno_order=5):
    """
    This example solves the 1-dimensional acoustics equations in a homogeneous
    medium.
    """
    import numpy as np

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
        solver.weno_order=weno_order
    else: raise Exception('Unrecognized value of solver_type.')

    #========================================================================
    # Instantiate the solver and define the system of equations to be solved
    #========================================================================
    solver.kernel_language=kernel_language
    from riemann import rp_acoustics
    solver.num_waves=rp_acoustics.num_waves
    if kernel_language=='Python': 
        solver.rp = rp_acoustics.rp_acoustics_1d
 
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.bc_lower[0] = pyclaw.BC.reflecting
    solver.bc_upper[0] = pyclaw.BC.reflecting

    solver.cfl_desired = 1.0
    solver.cfl_max     = 1.0

    #========================================================================
    # Instantiate the grid and set the boundary conditions
    #========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,200)
    grid = pyclaw.Grid(x)
    num_eqn = 2
    state = pyclaw.State(grid,num_eqn)

    #========================================================================
    # Set problem-specific variables
    #========================================================================
    rho = 1.0
    bulk = 1.0
    state.aux_global['rho']=rho
    state.aux_global['bulk']=bulk
    state.aux_global['zz']=np.sqrt(rho*bulk)
    state.aux_global['cc']=np.sqrt(bulk/rho)

    #========================================================================
    # Set the initial condition
    #========================================================================
    xc=grid.x.center
    state.q[0,:] = np.cos(2*np.pi*xc)
    state.q[1,:] = 0.
    
    #========================================================================
    # Set up the controller object
    #========================================================================
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir
    claw.nout = 40
    claw.tfinal = 2.0

    # Solve
    status = claw.run()

    # Plot results
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)
