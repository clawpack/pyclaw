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
    solver.mwaves = 2
    solver.mthlim = pyclaw.limiters.MC
    solver.kernel_language=kernel_language
    if kernel_language =='Python': 
        solver.set_riemann_solver('shallow_roe')
        grid.aux_global['g'] = 1.0
        grid.aux_global['efix'] = False
    solver.mthbc_lower[0] = pyclaw.BC.outflow
    solver.mthbc_upper[0] = pyclaw.BC.outflow

    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================
    xlower = -5.0
    xupper = 5.0
    mx = 500
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    grid = pyclaw.Grid(x)
    grid.meqn = 2

    # Parameters
    grid.aux_global['grav'] = 1.0

    grid.zeros_q()
    xCenter = grid.x.center

    damRadius = 0.0
    hl = 3.
    ul = 0.
    hr = 1.
    ur = 0.

    grid.q[0,:] = hl * (grid.p_center[0] <= damRadius) + hr * (grid.p_center[0] > damRadius)
    grid.q[1,:] = hl*ul * (grid.p_center[0] <= damRadius) + hr*ur * (grid.p_center[0] > damRadius)

    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = 2.0
    claw.solution = pyclaw.Solution(grid)
    claw.solver = solver
    claw.outdir = outdir


    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    #===========================================================================
    # Plot results
    #===========================================================================
    if iplot:     pyclaw.plot.plotInteractive(outdir=outdir,format=claw.output_format)
    if htmlplot:  pyclaw.plot.plotHTML(outdir=outdir,format=claw.output_format)


if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=shallow1D(*args,**kwargs)
    print 'Error: ',error
    

   

