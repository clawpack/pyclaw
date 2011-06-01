#!/usr/bin/env python
# encoding: utf-8

def burgers(use_petsc=False,kernel_language='Fortran',iplot=False,htmlplot=True,outdir='./_output',use_PETSc=False,solver_type='classic'):
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
    solver.mthlim = pyclaw.limiters.vanleer
    solver.mthbc_lower[0] = 2
    solver.mthbc_upper[0] = 2

    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,500)
    grid = pyclaw.Grid(x)
    grid.meqn = 1

    grid.zeros_q()
    xc=grid.x.center
    grid.q[0,:] = np.sin(np.pi*2*xc) + 0.50
    grid.aux_global['efix']=True

    #===========================================================================
    # Setup controller and controller parameters. Then solve the problem
    #===========================================================================
    output_format = 'petsc'
    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal =0.5
    claw.solution = pyclaw.Solution(grid)
    claw.solver = solver
    claw.outdir = outdir

    status = claw.run()

    if htmlplot:  pyclaw.plot.plotHTML(outdir=outdir)
    if iplot:     pyclaw.plot.plotInteractive(outdir=outdir)


if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=burgers(*args,**kwargs)
    print 'Error: ',error

