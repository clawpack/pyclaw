#!/usr/bin/env python
# encoding: utf-8
def advection(kernelsType='P',iplot=True,petscPlot=False,useController=True,usepetclaw=True,soltype='clawpack'):
    """
    Example python script for solving the 1d advection equation.
    """

    import numpy as np
    if usepetclaw:
        from petsc4py import PETSc
        from petclaw.grid import Grid, Dimension
        output_format='petsc'
        if soltype=='sharpclaw':
            from petclaw.evolve.sharpclaw import SharpPetClawSolver1D as ClawSolver
        else:
            from petclaw.evolve.clawpack import PetClawSolver1D as ClawSolver
        solver = ClawSolver(kernelsType = kernelsType)
    else: #Pure pyclaw
        output_format='ascii'
        from pyclaw.solution import Grid, Dimension
        if soltype=='sharpclaw':
            from pyclaw.evolve.sharpclaw import SharpClawSolver1D as ClawSolver
            solver = ClawSolver(kernelsType = kernelsType)
        else:
            from pyclaw.evolve.clawpack import ClawSolver1D as ClawSolver
            solver=ClawSolver()
        if kernelsType=='F': raise notImplementedError
    if soltype=='sharpclaw':
        solver.time_integrator='SSP33'
        solver.lim_type=2
        mbc=3
    else:
        mbc=2


    from pyclaw.solution import Solution
    from pyclaw.controller import Controller

    x = Dimension('x',0.0,1.0,200,mthbc_lower=2,mthbc_upper=2)
    x.mbc=mbc
    grid = Grid(x)
    grid.aux_global['u']=-1.
    if kernelsType=='F':
        from step1 import comrp
        comrp.u = grid.aux_global['u']
    grid.meqn = 1
    grid.t = 0.0
    grid.mbc=mbc

    if usepetclaw: grid.init_q_petsc_structures()

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
    solver.cfl_max=0.25
    solver.cfl_desired=0.15
    if kernelsType=='P': solver.set_riemann_solver('advection')

    claw.keep_copy = False
    claw.dt_variable = True
    claw.outstyle=1
    claw.nout=10
    claw.tfinal =2.0
    claw.solutions['n'] = initial_solution
    claw.solver = solver

    status = claw.run()

    from petclaw import plot
    plot.plotInteractive(format=output_format)

    output_object=claw

    print np.max(np.abs((claw.solutions['n'].grids[0].q-q0)))*grid.d[0]
    return output_object

if __name__=="__main__":
    advection()
