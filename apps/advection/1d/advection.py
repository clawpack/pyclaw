#!/usr/bin/env python
# encoding: utf-8
def advection(kernel_language='Python',iplot=True,petscPlot=False,useController=True,use_PETSc=False,soltype='classic'):
    """
    Example python script for solving the 1d advection equation.
    """

    import numpy as np
    if use_PETSc:
        from petsc4py import PETSc
        from petclaw.grid import Grid, Dimension
        output_format='petsc'
        if soltype=='sharpclaw':
            from petclaw.evolve.sharpclaw import SharpPetClawSolver1D as ClawSolver
        else:
            from petclaw.evolve.clawpack import PetClawSolver1D as ClawSolver
    else: #Pure pyclaw
        output_format='ascii'
        from pyclaw.grid import Grid, Dimension
        if soltype=='sharpclaw':
            from pyclaw.evolve.sharpclaw import SharpClawSolver1D as ClawSolver
        else:
            from pyclaw.evolve.clawpack import ClawSolver1D as ClawSolver
        solver = ClawSolver()
        if kernel_language=='Fortran': raise notImplementedError
    if soltype=='sharpclaw':
        solver.time_integrator='SSP33'
        solver.lim_type=2
        mbc=3
    else:
        mbc=2


    from pyclaw.solution import Solution
    from pyclaw.controller import Controller

    x = Dimension('x',0.0,1.0,50,mthbc_lower=2,mthbc_upper=2)
    x.mbc=mbc
    grid = Grid(x)
    grid.aux_global['u']=-1.
    grid.meqn = 1
    grid.t = 0.0
    grid.mbc=mbc

    if use_PETSc: grid.init_q_petsc_structures()

    xc=grid.x.center
    beta=100; gamma=0; x0=0.75
    q=np.zeros([grid.meqn,len(xc)],order='F')
    q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    grid.q = q
    q0 = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

    initial_solution = Solution(grid)

    claw = Controller()

    claw.output_format=output_format
    solver.mwaves = 1
    solver.cfl_max=0.5
    solver.cfl_desired=0.45
    solver.kernel_language='Python'
    if solver.kernel_language=='Python': 
        solver.set_riemann_solver('advection')
    else:
        from step1 import comrp
        comrp.u = grid.aux_global['u']

    claw.keep_copy = False
    claw.dt_variable = True
    claw.outstyle=1
    claw.nout=10
    claw.tfinal =2.0
    claw.solutions['n'] = initial_solution
    claw.solver = solver

    status = claw.run()

    if iplot:
        from petclaw import plot
        plot.plotInteractive(format=output_format)

    output_object=claw

    print np.max(np.abs((claw.solutions['n'].grids[0].q-q0)))*grid.d[0]
    return output_object

if __name__=="__main__":
    import sys
    from pyclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    output=advection(*args,**kwargs)
