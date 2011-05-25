#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(use_PETSc=False,kernel_language='Fortran',soltype='classic',iplot=False,htmlplot=False,outdir='./_output'):
    """
    This example solves the 1-dimensional acoustics equations in a homogeneous
    medium.
    """
    import numpy as np
    from petclaw import plot

    #======================================================
    # Import the appropriate classes
    #======================================================
    if use_PETSc:
        import petclaw as myclaw
        if soltype=='classic':
            from petclaw.evolve.clawpack import PetClawSolver1D as mySolver
        elif soltype=='sharpclaw':
            from petclaw.evolve.sharpclaw import PetSharpClawSolver1D as mySolver
        else: raise Exception('Unrecognized value of soltype.')
        from petclaw.controller import Controller
    else: #Pure pyclaw
        import pyclaw as myclaw
        if soltype=='classic':
            from pyclaw.evolve.clawpack import ClawSolver1D as mySolver
        elif soltype=='sharpclaw':
            from pyclaw.evolve.sharpclaw import SharpClawSolver1D as mySolver
        else: raise Exception('Unrecognized value of soltype.')
        from pyclaw.controller import Controller

    from pyclaw.solution import Solution

    #========================================================================
    # Instantiate the solver and define the system of equations to be solved
    #========================================================================
    solver = mySolver()
    solver.mwaves=2
    solver.kernel_language=kernel_language
    if kernel_language=='Python': solver.set_riemann_solver('acoustics')
 
    from pyclaw.evolve import limiters
    solver.mthlim = limiters.MC

    #========================================================================
    # Instantiate the grid
    #========================================================================
    x = myclaw.grid.Dimension('x',0.0,1.0,100,mthbc_lower=3,mthbc_upper=1)
    grid = myclaw.grid.Grid(x)
    grid.meqn=2
    grid.mbc=solver.mbc
    # init_q_petsc_structures must be called 
    # before grid.x.center and such can be accessed.
    if use_PETSc: grid.init_q_petsc_structures()

    #========================================================================
    # This part should really just depend on the solver, not the grid
    #========================================================================
    rho = 1.0
    bulk = 1.0
    grid.aux_global['rho']=rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']=np.sqrt(rho*bulk)
    grid.aux_global['cc']=np.sqrt(rho/bulk)
    if kernel_language=='Fortran':
        if soltype=='classic':
            from classic1 import cparam 
        elif soltype=='sharpclaw':
            from sharpclaw1 import cparam
        for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)


    #========================================================================
    # Set the initial condition
    #========================================================================
    xc=grid.x.center
    q=np.zeros([grid.meqn,len(xc)], order = 'F')
    beta=100; gamma=0; x0=0.75
    q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    q[1,:]=0.
    grid.q=q
    
    #========================================================================
    # Set up the controller object
    #========================================================================
    claw = Controller()

    claw.solutions['n'] = Solution(grid)
    claw.solver = solver

    claw.nout = 10
    claw.outdir = outdir
    claw.tfinal = 1.0

    # Solve
    status = claw.run()

    # Plot results
    if htmlplot:  plot.plotHTML(outdir=outdir,format=claw.output_format)
    if iplot:     plot.plotInteractive(outdir=outdir,format=claw.output_format)

if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=acoustics(*args,**kwargs)
