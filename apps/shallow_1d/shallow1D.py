#!/usr/bin/env python
# encoding: utf-8

"""
1D shallow water equations.
"""

    
def shallow1D(use_petsc=False,kernel_language='Fortran',iplot=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D()
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D()

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    solver.num_waves = 2
    solver.limiters = pyclaw.limiters.tvd.vanleer
    solver.kernel_language=kernel_language
    if kernel_language =='Python': 
        solver.set_riemann_solver('shallow_roe')
        grid.problem_data['g'] = 1.0
        grid.problem_data['efix'] = False
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================
    xlower = -5.0
    xupper = 5.0
    mx = 500
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    grid = pyclaw.Grid(x)
    num_eqn = 2
    state = pyclaw.State(grid,num_eqn)

    # Parameters
    state.problem_data['grav'] = 1.0

    xc = grid.x.center

    IC='2-shock'
    x0=0.

    if IC=='dam-break':
        hl = 3.
        ul = 0.
        hr = 1.
        ur = 0.
        state.q[0,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[1,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='2-shock':
        hl = 1.
        ul = 1.
        hr = 1.
        ur = -1.
        state.q[0,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[1,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='perturbation':
        eps=0.1
        state.q[0,:] = 1.0 + eps*np.exp(-(xc-x0)**2/0.5)
        state.q[1,:] = 0.

    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = True
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
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir,format=claw.output_format)
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir,format=claw.output_format)


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow1D)
    

   

