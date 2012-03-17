#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(use_petsc=True,kernel_language='Fortran',solver_type='classic',iplot=False,htmlplot=False,outdir='./_output',weno_order=5):
    import numpy as np
    """
    1D acoustics example.
    """

    output_format=None #Suppress output to make tests faster
    if use_petsc:
        import petclaw as pyclaw
    else: #Pure pyclaw
        import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
        solver.weno_order = weno_order
    else: raise Exception('Unrecognized value of solver_type.')

    # Initialize patches and solution
    x = pyclaw.Dimension('x',0.0,1.0,100)
    domain = pyclaw.Domain(x)
    num_eqn=2
    state = pyclaw.State(domain,num_eqn)

    rho = 1.0
    bulk = 1.0
    state.problem_data['rho']=rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']=np.sqrt(rho*bulk)
    state.problem_data['cc']=np.sqrt(rho/bulk)

    grid = state.grid
    xc=grid.x.centers
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    state.q[1,:]=0.
    
    init_solution = pyclaw.Solution(state,domain)


    solver.num_waves=2
    solver.kernel_language=kernel_language

    if kernel_language=='Python': 
        from riemann import rp_acoustics
        solver.rp = rp_acoustics.rp_acoustics_1d
    elif kernel_language=='Fortran':
        import riemann
        solver.rp = riemann.rp1_acoustics


    solver.limiters = [4]*solver.num_waves
    solver.dt_initial=grid.delta[0]/state.problem_data['cc']*0.1
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.num_output_times = 5
    
    claw.output_format = output_format

    claw.outdir = outdir
    claw.tfinal = 1.0
    claw.solution = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    #This test is set up so that the waves pass through the domain
    #exactly once, and the final solution should be equal to the
    #initial condition.  Here we output the 1-norm of their difference.
    if use_petsc==True:
        q0=claw.frames[0].state.gqVec.getArray().reshape([-1])
        qfinal=claw.frames[claw.num_output_times].state.gqVec.getArray().reshape([-1])
    else:
        q0=claw.frames[0].state.q.reshape([-1])
        qfinal=claw.frames[claw.num_output_times].state.q.reshape([-1])
    dx=claw.solution.domain.grid.delta[0]

    return dx*np.sum(np.abs(qfinal-q0))


if __name__=="__main__":
    import sys
    from pyclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=acoustics(*args,**kwargs)
    print '1-norm of difference between initial and final solution: ',error
