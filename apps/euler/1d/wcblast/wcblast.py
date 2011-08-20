#!/usr/bin/env python
# encoding: utf-8
r"""Woodward-Colella blast wave interaction problem.
   1D compressible inviscid flow (Euler equations)."""

gamma = 1.4
gamma1 = gamma - 1.

def wcblast(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a pair of interacting shock waves.
    The conserved quantities are density, momentum density, and total energy density.
    """

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()

    solver.mwaves = 3
    solver.mthbc_lower[0]=pyclaw.BC.reflecting
    solver.mthbc_upper[0]=pyclaw.BC.reflecting

    # Initialize grid
    mx=500;
    x = pyclaw.Dimension('x',0.0,1.0,mx)
    grid = pyclaw.Grid([x])
    meqn = 3
    state = pyclaw.State(grid,meqn)

    state.aux_global['gamma']= gamma
    state.aux_global['gamma1']= gamma1

    state.q[0,:] = 1.
    state.q[1,:] = 0.
    x =state.grid.x.center
    state.q[2,:] = ( (x<0.1)*1.e3 + (0.1<=x)*(x<0.9)*1.e-2 + (0.9<=x)*1.e2 ) / gamma1

    solver.limiters = 4

    claw = pyclaw.Controller()
    claw.tfinal = 0.038
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.nout = 10
    claw.outdir = outdir

    # Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    return claw.solution.q

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(wcblast)
