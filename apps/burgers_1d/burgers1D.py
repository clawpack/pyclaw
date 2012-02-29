#!/usr/bin/env python
# encoding: utf-8

def burgers(use_petsc=0,kernel_language='Fortran',iplot=0,htmlplot=0,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 1d Burgers equation.
    """

    import numpy as np

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()
        solver.limiters = pyclaw.limiters.tvd.vanleer

    solver.kernel_language = kernel_language
    if kernel_language=='Python': 
        solver.set_riemann_solver('burgers')
    elif kernel_language=='Fortran':
        import riemann
        solver.rp = riemann.rp1_burgers
        
    solver.num_waves = 1
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    #===========================================================================
    # Initialize domain and then initialize the solution associated to the domain
    #===========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,500)
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    grid = state.grid
    xc=grid.x.centers
    state.q[0,:] = np.sin(np.pi*2*xc) + 0.50
    state.problem_data['efix']=True

    #===========================================================================
    # Setup controller and controller parameters. Then solve the problem
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal =0.5
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir

    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(burgers)

