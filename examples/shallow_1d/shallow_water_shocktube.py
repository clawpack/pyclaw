#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow
==================

Solve the one-dimensional shallow water equations:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = 0.

Here h is the depth, u is the velocity, and g is the gravitational constant.
"""

    
def setup(use_petsc=False,kernel_language='Fortran',outdir='./_output',solver_type='classic'):
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language =='Python':
        rs = riemann.shallow_1D_py.shallow_1D
    elif kernel_language =='Fortran':
        rs = riemann.shallow_roe_with_efix_1D
 
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    solver.kernel_language=kernel_language

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    #===========================================================================
    # Initialize domain and then initialize the solution associated to the domain
    #===========================================================================
    xlower = -5.0
    xupper = 5.0
    mx = 500
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)

    # Parameters
    state.problem_data['grav'] = 1.0
    
    xc = state.grid.x.centers

    IC='2-shock'
    x0=0.

    if IC=='dam-break':
        hl = 3.
        ul = 0.
        hr = 1.
        ur = 0.
        state.q[0,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[1,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='2-shock':
        hl = 1.
        ul = 1.
        hr = 1.
        ur = -1.
        state.q[0,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[1,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='perturbation':
        eps=0.1
        state.q[0,:] = 1.0 + eps*np.exp(-(xc-x0)**2/0.5)
        state.q[1,:] = 0.

    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = 2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir

    return claw


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
