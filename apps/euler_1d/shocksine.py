#!/usr/bin/env python
# encoding: utf-8
r"""Shu-Osher problem.
   1D compressible inviscid flow (Euler equations)."""

gamma = 1.4
gamma1 = gamma - 1.

def soshock(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a shock wave impacting a sinusoidal density field.
    """
    import numpy as np

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()
        solver.limiters = 4

    from clawpack import riemann
    solver.rp = riemann.rp1_euler_with_efix

    solver.num_waves = 3
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap

    # Initialize domain
    mx=400;
    x = pyclaw.Dimension('x',-5.0,5.0,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 3
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1

    xc =state.grid.x.centers
    epsilon=0.2
    state.q[0,:] = (xc<-4.)*3.857143 + (xc>=-4.)*(1+epsilon*np.sin(5*xc))
    velocity = (xc<-4.)*2.629369
    state.q[1,:] = velocity * state.q[0,:]
    pressure = (xc<-4.)*10.33333 + (xc>=-4.)*1.
    state.q[2,:] = pressure/gamma1 + 0.5 * state.q[0,:] * velocity**2

    claw = pyclaw.Controller()
    claw.tfinal = 1.8
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.outdir = outdir

    # Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    return claw.solution.q

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(soshock)
