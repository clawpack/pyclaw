#!/usr/bin/env python
# encoding: utf-8
"""
Example python script for solving the 1d variable-coefficient advection 
equation: q_t + u(x)q_x = 0.
"""

import numpy as np

def qinit(state):

    # Initial Data parameters
    ic = 3
    beta = 100.
    gamma = 0.
    x0 = 0.3
    x1 = 0.7
    x2 = 0.9

    x =state.grid.x.centers
    
    # Gaussian
    qg = np.exp(-beta * (x-x0)**2) * np.cos(gamma * (x - x0))
    # Step Function
    qs = (x > x1) * 1.0 - (x > x2) * 1.0
    
    if   ic == 1: state.q[0,:] = qg
    elif ic == 2: state.q[0,:] = qs
    elif ic == 3: state.q[0,:] = qg + qs


def auxinit(state):
    # Initilize petsc Structures for aux
    xc=state.grid.x.centers
    state.aux[0,:] = np.sin(2.*np.pi*xc)+2
    

def setup(use_petsc=False,solver_type='classic',kernel_language='Python',outdir='./_output'):
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        if kernel_language == 'Fortran':
            solver = pyclaw.ClawSolver1D(riemann.vc_advection_1D)
        elif kernel_language=='Python': 
            solver = pyclaw.ClawSolver1D(riemann.vc_advection_1D_py.vc_advection_1D)
    elif solver_type=='sharpclaw':
        if kernel_language == 'Fortran':
            solver = pyclaw.SharpClawSolver1D(riemann.vc_advection_1D)
        elif kernel_language=='Python': 
            solver = pyclaw.SharpClawSolver1D(riemann.vc_advection_1D_py.vc_advection_1D)
        solver.weno_order=weno_order
    else: raise Exception('Unrecognized value of solver_type.')

    solver.kernel_language = kernel_language

    solver.limiters = pyclaw.limiters.tvd.MC
    solver.bc_lower[0] = 2
    solver.bc_upper[0] = 2
    solver.aux_bc_lower[0] = 2
    solver.aux_bc_upper[0] = 2

    xlower=0.0; xupper=1.0; mx=100
    x    = pyclaw.Dimension('x',xlower,xupper,mx)
    domain = pyclaw.Domain(x)
    num_aux=1
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn,num_aux)

    qinit(state)
    auxinit(state)

    claw = pyclaw.Controller()
    claw.outdir = outdir
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    claw.tfinal = 1.0
    
    return claw


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
