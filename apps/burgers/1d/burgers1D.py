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

    solver.kernel_language = kernel_language
    if kernel_language=='Python': solver.set_riemann_solver('burgers')
    solver.mwaves = 1
    solver.limiters = pyclaw.limiters.tvd.vanleer
    solver.mthbc_lower[0] = pyclaw.BC.periodic
    solver.mthbc_upper[0] = pyclaw.BC.periodic

    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,500)
    grid = pyclaw.Grid(x)
    meqn = 1
    state = pyclaw.State(grid,meqn)

    xc=grid.x.center
    state.q[0,:] = np.sin(np.pi*2*xc) + 0.50
    state.aux_global['efix']=True

    #===========================================================================
    # Setup controller and controller parameters. Then solve the problem
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal =0.5
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir

    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(burgers)
    print 'Error: ',output

