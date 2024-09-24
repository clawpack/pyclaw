#!/usr/bin/env python
# encoding: utf-8

def acoustics(problem='figure 9.4'):
    """
    This example solves the 1-dimensional variable-coefficient acoustics
    equations in a medium with a single interface.
    """
    from numpy import sqrt, abs
    from clawpack import pyclaw
    from clawpack import riemann

    solver = pyclaw.ClawSolver1D(riemann.acoustics_variable_1D)

    solver.limiters = pyclaw.limiters.tvd.MC
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

    x = pyclaw.Dimension(-5.0,5.0,500,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 2
    num_aux = 2
    state = pyclaw.State(domain, num_eqn, num_aux)

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
    xc = domain.grid.x.centers

    state.aux[0,:] = (xc<=0)*zl + (xc>0)*zr  # Impedance
    state.aux[1,:] = (xc<=0)*cl + (xc>0)*cr  # Sound speed

    # initial condition: half-ellipse
    state.q[0,:] = sqrt(abs(1.-(xc+3.)**2))*(xc>-4.)*(xc<-2.)
    state.q[1,:] = state.q[0,:] + 0.

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.tfinal = 5.0
    claw.num_output_times   = 10

    # Solve
    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)
