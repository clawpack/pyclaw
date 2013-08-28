#!/usr/bin/env python
# encoding: utf-8
r"""Shu-Osher problem.
   1D compressible inviscid flow (Euler equations)."""

import numpy as np
gamma = 1.4
gamma1 = gamma - 1.

a = np.array([[0., 0., 0., 0., 0., 0., 0.],
              [.3772689153313680, 0., 0., 0., 0., 0., 0.],
              [.3772689153313680, .3772689153313680, 0., 0., 0., 0., 0.],
              [.2429952205373960, .2429952205373960, .2429952205373960, 0., 0., 0., 0.],
              [.1535890676951260, .1535890676951260, .1535890676951260, .2384589328462900, 0., 0., 0.]])

c = np.array([0., .3772689153313680, .7545378306627360, .7289856616121880, .6992261359316680])

b = np.array([.206734020864804, .206734020864804, .117097251841844, .181802560120140, .287632146308408])

def soshock(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='sharpclaw'):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a shock wave impacting a sinusoidal density field.
    """
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann.euler_with_efix_1D)
        solver.time_integrator = 'RK'
        solver.a, solver.b, solver.c = a, b, c
        solver.cfl_desired = 0.6
        solver.cfl_max = 0.7
    else:
        solver = pyclaw.ClawSolver1D(riemann.euler_with_efix_1D)

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap

    # Initialize domain
    mx=400;
    x = pyclaw.Dimension('x',-5.0,5.0,mx)
    domain = pyclaw.Domain([x])
    state = pyclaw.State(domain,solver.num_eqn)

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
