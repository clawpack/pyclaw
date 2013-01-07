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
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D()
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D()

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    solver.kernel_language=kernel_language
    from clawpack.riemann import rp_shallow
    solver.num_waves = rp_shallow.num_waves
    if kernel_language =='Python':
        solver.rp = rp_shallow.rp_shallow_roe_1d
        state.problem_data['efix'] = True
    elif kernel_language == 'Fortran':
        from clawpack import riemann
        solver.rp = riemann.rp1_shallow_roe_with_efix

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    #===========================================================================
    # Initialize domain and then initialize the solution associated to the domain
    #===========================================================================
    xlower = -5.0
    xupper = 5.0
    mx = 500
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)

    # Parameters
    state.problem_data['grav'] = 1.0
    
    xc = state.grid.x.centers

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
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir


    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    #===========================================================================
    # Plot results
    #===========================================================================
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir,file_format=claw.output_format)
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir,file_format=claw.output_format)


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow1D)
    

   

