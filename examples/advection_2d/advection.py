#!/usr/bin/env python
# encoding: utf-8

#===========================================================================
# Import libraries
#===========================================================================
import numpy as np

def qinit(state):

    # Set initial conditions for q.
    # Sample scalar equation with data that is piecewise constant with
    # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    #     0.1  otherwise
    
    x = state.grid.x.centers
    y = state.grid.y.centers
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] > 0.0 and x[i] < 0.5 and y[j]>0.0 and y[j] < 0.5:
                state.q[:,i,j] = 1.0
            else:
                state.q[:,i,j] = 0.1
                
def advection2D(use_petsc=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d advection equation.
    """
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D(riemann.advection_2D)
        solver.dimensional_split = 1
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.advection_2D)

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic

    solver.cfl_max=1.0
    solver.cfl_desired = 0.9

    #===========================================================================
    # Initialize domain, then initialize the solution associated to the domain and
    # finally initialize aux array
    #===========================================================================

    # Domain:
    mx=50; my=50
    x = pyclaw.Dimension('x',0.0,1.0,mx)
    y = pyclaw.Dimension('y',0.0,1.0,my)
    domain = pyclaw.Domain([x,y])

    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['u'] = 0.5 # Parameters (global auxiliary variables)
    state.problem_data['v'] = 1.0

    # Initial solution
    # ================
    qinit(state) # This function is defined above


    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir

    return claw


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(advection2D)
