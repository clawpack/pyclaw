#!/usr/bin/env python
# encoding: utf-8
r"""Woodward-Colella blast wave interaction problem.
   1D compressible inviscid flow (Euler equations)."""

gamma = 1.4
gamma1 = gamma - 1.

def setup(use_petsc=False,outdir='./_output',solver_type='classic'):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a pair of interacting shock waves.
    The conserved quantities are density, momentum density, and total energy density.
    """
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann.euler_with_efix_1D)
    else:
        solver = pyclaw.ClawSolver1D(riemann.euler_with_efix_1D)

    solver.bc_lower[0]=pyclaw.BC.wall
    solver.bc_upper[0]=pyclaw.BC.wall

    # Initialize domain
    mx=500;
    x = pyclaw.Dimension('x',0.0,1.0,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 3
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1

    state.q[0,:] = 1.
    state.q[1,:] = 0.
    x =state.grid.x.centers
    state.q[2,:] = ( (x<0.1)*1.e3 + (0.1<=x)*(x<0.9)*1.e-2 + (0.9<=x)*1.e2 ) / gamma1

    solver.limiters = 4

    claw = pyclaw.Controller()
    claw.tfinal = 0.038
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.outdir = outdir

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
