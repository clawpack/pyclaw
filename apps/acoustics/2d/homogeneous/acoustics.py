#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,width=0.2):
    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    q[1,:,:] = 0.
    q[2,:,:] = 0.
    grid.q=q


def acoustics2D(iplot=False,htmlplot=False,use_PETSc=False,outdir='./_output',soltype='classic'):
    """
    Example python script for solving the 2d acoustics equations.
    """

    if use_PETSc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if soltype=='classic':
        solver=pyclaw.evolve.clawpack.ClawSolver2D()
    elif soltype=='sharpclaw':
        solver=pyclaw.evolve.clawpack.SharpClawSolver2D()

    from petclaw import plot

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 2
    from pyclaw.evolve import limiters
    solver.mthlim = limiters.MC
    solver.dim_split=False

    solver.mthbc_lower[0]=pyclaw.BC.outflow
    solver.mthbc_upper[0]=pyclaw.BC.outflow
    solver.mthbc_lower[1]=pyclaw.BC.outflow
    solver.mthbc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=100; my=100
    x = pyclaw.Dimension('x',-1.0,1.0,mx,mthbc_lower=1,mthbc_upper=1)
    y = pyclaw.Dimension('y',-1.0,1.0,my,mthbc_lower=1,mthbc_upper=1)
    grid = pyclaw.Grid([x,y])

    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    grid.aux_global['rho']= rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']= zz
    grid.aux_global['cc']=cc

    solver.dt=np.min(grid.d)/grid.aux_global['cc']*solver.cfl_desired

    grid.meqn = 3
    grid.mbc = solver.mbc
    tfinal = 0.12

    if use_PETSc:
        # Initialize petsc Structures for q
        grid.init_q_petsc_structures()

    qinit(grid)
    initial_solution = pyclaw.Solution(grid)

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = tfinal
    claw.solutions['n'] = initial_solution
    claw.solver = solver
    claw.outdir=outdir

    # Solve
    status = claw.run()

    if htmlplot:  plot.plotHTML(outdir=outdir,format=claw.output_format)
    if iplot:     plot.plotInteractive(outdir=outdir,format=claw.output_format)

    if use_PETSc:
        pressure=claw.frames[claw.nout].grid.gqVec.getArray().reshape([grid.local_n[0],grid.local_n[1],grid.meqn])[:,:,0]
    else:
        pressure=claw.frames[claw.nout].grid.q[:,:,0]
    return pressure


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        acoustics2D(*args,**kwargs)
    else: acoustics2D()
