#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(use_petsc=True,kernel_language='Fortran',solver_type='classic',iplot=False,htmlplot=False,outdir='./_output'):
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
    else: raise Exception('Unrecognized value of solver_type.')

    # Initialize grids and solution
    x = pyclaw.Dimension('x',0.0,1.0,100)
    grid = pyclaw.Grid(x)
    meqn=2
    state = pyclaw.State(grid,meqn)

    rho = 1.0
    bulk = 1.0
    state.aux_global['rho']=rho
    state.aux_global['bulk']=bulk
    state.aux_global['zz']=np.sqrt(rho*bulk)
    state.aux_global['cc']=np.sqrt(rho/bulk)

    xc=grid.x.center
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    state.q[1,:]=0.
    
    init_solution = pyclaw.Solution(state)


    solver.mwaves=2
    solver.kernel_language=kernel_language

    if kernel_language=='Python': 
        from riemann import rp_acoustics
        solver.rp = rp_acoustics.rp_acoustics_1d

    solver.limiters = [4]*solver.mwaves
    solver.dt=grid.d[0]/state.aux_global['cc']*0.1
    solver.mthbc_lower[0] = pyclaw.BC.periodic
    solver.mthbc_upper[0] = pyclaw.BC.periodic

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.nout = 5
    
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
        qfinal=claw.frames[claw.nout].state.gqVec.getArray().reshape([-1])
    else:
        q0=claw.frames[0].state.q.reshape([-1])
        qfinal=claw.frames[claw.nout].state.q.reshape([-1])
    dx=claw.frames[0].grid.d[0]

    return dx*np.sum(np.abs(qfinal-q0))


if __name__=="__main__":
    import sys
    from pyclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=acoustics(*args,**kwargs)
    print '1-norm of difference between initial and final solution: ',error
