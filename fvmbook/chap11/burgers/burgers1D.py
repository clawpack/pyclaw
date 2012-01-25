#!/usr/bin/env python
# encoding: utf-8

def burgers(iplot=1,htmlplot=0,outdir='./_output'):
    """
    Example from Chapter 11 of LeVeque, Figure 11.8.
    Shows decay of an initial wave packet to an N-wave with Burgers' equation.
    """
    import numpy as np

    import pyclaw

    solver = pyclaw.ClawSolver1D()

    solver.num_waves = 1
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================
    x = pyclaw.Dimension('x',-8.0,8.0,1000)
    grid = pyclaw.Grid(x)
    num_eqn = 1
    state = pyclaw.State(grid,num_eqn)

    xc=grid.x.center
    state.q[0,:] = (xc>-np.pi)*(xc<np.pi)*(2.*np.sin(3.*xc)+np.cos(2.*xc)+0.2)
    state.q[0,:] = state.q[0,:]*(np.cos(xc)+1.)
    state.aux_global['efix']=True

    #===========================================================================
    # Setup controller and controller parameters. Then solve the problem
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 6.0
    claw.nout   = 30
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir

    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(burgers)
