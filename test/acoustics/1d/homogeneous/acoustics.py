#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(use_PETSc=True,kernel_language='Fortran',soltype='classic',iplot=False,htmlplot=False,outdir='./_output'):
    import numpy as np
    """
    1D acoustics example.
    """

    output_format=None #Suppress output to make tests faster
    if use_PETSc:
        import petclaw as pyclaw
    else: #Pure pyclaw
        import pyclaw

    if soltype=='classic':
        from pyclaw.evolve.clawpack import ClawSolver1D 
        solver = ClawSolver1D()
    elif soltype=='sharpclaw':
        from pyclaw.evolve.sharpclaw import SharpClawSolver1D
        solver = SharpClawSolver1D()
    else: raise Exception('Unrecognized value of soltype.')

    from pyclaw.solution import Solution
    from pyclaw.controller import Controller

    # Initialize grids and solutions
    x = pyclaw.grid.Dimension('x',0.0,1.0,100,mthbc_lower=2,mthbc_upper=2)
    grid = pyclaw.grid.Grid(x)
    grid.meqn=2

    if soltype=='classic': grid.mbc=2
    elif soltype=='sharpclaw': grid.mbc=3

    rho = 1.0
    bulk = 1.0
    grid.aux_global['rho']=rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']=np.sqrt(rho*bulk)
    grid.aux_global['cc']=np.sqrt(rho/bulk)

    # init_q_petsc_structures must be called 
    # before grid.x.center and such can be accessed.
    if use_PETSc:
        grid.init_q_petsc_structures()

    xc=grid.x.center
    q=np.zeros([grid.meqn,len(xc)], order = 'F')
    beta=100; gamma=0; x0=0.75
    q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    q[1,:]=0.
    grid.q=q
    
    init_solution = Solution(grid)


    solver.mwaves=2
    solver.kernel_language=kernel_language

    if kernel_language=='Python': 
        solver.set_riemann_solver('acoustics')

    solver.mthlim = [4]*solver.mwaves
    solver.dt=grid.d[0]/grid.aux_global['cc']*0.1

    claw = Controller()
    claw.keep_copy = True
    claw.nout = 5
    
    claw.output_format = output_format

    claw.outdir = outdir
    claw.tfinal = 1.0
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    #This test is set up so that the waves pass through the domain
    #exactly once, and the final solution should be equal to the
    #initial condition.  Here we output the 1-norm of their difference.
    if use_PETSc==True:
        q0=claw.frames[0].grid.gqVec.getArray().reshape([-1])
        qfinal=claw.frames[claw.nout].grid.gqVec.getArray().reshape([-1])
    else:
        q0=claw.frames[0].grid.q.reshape([-1])
        qfinal=claw.frames[claw.nout].grid.q.reshape([-1])
    dx=claw.frames[0].grid.d[0]

    return dx*np.sum(np.abs(qfinal-q0))


if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=acoustics(*args,**kwargs)
    print '1-norm of difference between initial and final solutions: ',error
