#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(kernelsType='F',petscPlot=False,iplot=False,htmlplot=False,outdir='./_output'):
    import numpy as np
    """
    1D acoustics example.
    """

    from pyclaw.solution import Grid
    from pyclaw.solution import Dimension
    from pyclaw.solution import Solution
    from pyclaw.evolve.sharpclaw import SharpClawSolver1D
    from pyclaw.controller import Controller
    from petclaw import plot


    # Initialize grids and solutions
    x = Dimension('x',0.0,1.0,100,mthbc_lower=2,mthbc_upper=2)
    grid = Grid(x)
    rho = 1.0
    bulk = 1.0
    grid.aux_global['rho']=rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']=np.sqrt(rho*bulk)
    grid.aux_global['cc']=np.sqrt(rho/bulk)
    from flux1 import cparam 
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)
    grid.meqn=2
    grid.mbc=3
    grid.t = 0.0

    # init_q_petsc_structures must be called 
    # before grid.x.center and such can be accessed.
    xc=grid.x.center
    q=np.zeros([grid.meqn,len(xc)], order = 'F')
    beta=100; gamma=0; x0=0.75
    q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    q[1,:]=0.
    grid.q=q
    
    init_solution = Solution(grid)

    solver = SharpClawSolver1D(kernelsType = kernelsType)
    solver.mwaves=2
    if kernelsType=='P': solver.set_riemann_solver('acoustics')
    solver.mthlim = [4]*solver.mwaves
    solver.dt=grid.d[0]/grid.aux_global['cc']*0.1
    solver.time_integrator='SSP33'
    solver.cfl_desired=0.45
    solver.cfl_max=0.5

    claw = Controller()
    claw.keep_copy = True
    claw.outstyle = 3
    claw.nout = 5
    claw.iout = 5
    # The output format MUST be set to petsc!
    claw.output_format = 'ascii'
    claw.outdir = outdir
    claw.tfinal = 1.0
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    if htmlplot:  plot.plotHTML()
    if petscPlot: plot.plotPetsc(output_object)
    if iplot:     plot.plotInteractive(format=claw.output_format)

    #This test is set up so that the waves pass through the domain
    #exactly once, and the final solution should be equal to the
    #initial condition.  Here we output the 1-norm of their difference.
    #q0=claw.frames[0].grid.gqVec.getArray().reshape([-1])
    #qfinal=claw.frames[claw.nout].grid.gqVec.getArray().reshape([-1])
    #dx=claw.frames[0].grid.d[0]

    return 1#dx*np.sum(np.abs(qfinal-q0))


if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=acoustics(*args,**kwargs)
    print 'Error: ',error
