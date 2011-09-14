#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(solver_type='classic',iplot=False,htmlplot=False,outdir='./_output',problem='figure 9.4'):
    """
    This example solves the 1-dimensional variable-coefficient acoustics
    equations in a medium with a single interface.
    """
    from numpy import sqrt, abs

    import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else: raise Exception('Unrecognized value of solver_type.')

    solver.mwaves=2
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.mthbc_lower[0] = pyclaw.BC.outflow
    solver.mthbc_upper[0] = pyclaw.BC.outflow
    solver.mthauxbc_lower[0] = pyclaw.BC.outflow
    solver.mthauxbc_upper[0] = pyclaw.BC.outflow

    x = pyclaw.Dimension('x',-5.0,5.0,500)
    grid = pyclaw.Grid(x)
    meqn = 2
    maux = 2
    state = pyclaw.State(grid,meqn,maux)

    if problem == 'figure 9.4':
        rhol = 1.0
        cl   = 1.0
        rhor = 2.0
        cr   = 0.5
    elif problem == 'figure 9.5':
        rhol = 1.0
        cl   = 1.0
        rhor = 4.0
        cr   = 0.5
    zl = rhol*cl
    zr = rhor*cr
    xc = grid.x.center

    state.aux[0,:] = (xc<=0)*zl + (xc>0)*zr  # Impedance
    state.aux[1,:] = (xc<=0)*cl + (xc>0)*cr  # Sound speed

    # initial condition: half-ellipse
    state.q[0,:] = sqrt(abs(1.-(xc+3.)**2))*(xc>-4.)*(xc<-2.)
    state.q[1,:] = state.q[0,:] + 0.

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.tfinal = 5.0
    claw.nout   = 10

    # Solve
    status = claw.run()

    # Plot results
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)
