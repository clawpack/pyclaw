#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def acoustics2D(iplot=False,htmlplot=False,use_PETSc=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d acoustics equations.
    """

    if use_PETSc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver=pyclaw.evolve.clawpack.ClawSolver2D()
    elif solver_type=='sharpclaw':
        solver=pyclaw.evolve.clawpack.SharpClawSolver2D()

    solver.dim_split=False
    solver.mwaves = 2
    solver.mthlim = pyclaw.limiters.MC

    solver.mthbc_lower[0]=pyclaw.BC.reflecting
    solver.mthbc_upper[0]=pyclaw.BC.outflow
    solver.mthbc_lower[1]=pyclaw.BC.reflecting
    solver.mthbc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=100; my=100
    x = pyclaw.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.Dimension('y',-1.0,1.0,my)
    grid = pyclaw.Grid([x,y])

    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    grid.aux_global['rho']= rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']= zz
    grid.aux_global['cc']=cc

    grid.meqn = 3
    grid.mbc = solver.mbc

    Y,X = np.meshgrid(grid.y.center,grid.x.center)
    r = np.sqrt(X**2 + Y**2)
    width=0.2
    grid.zeros_q()
    grid.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(grid)
    claw.solver = solver
    claw.outdir=outdir

    # Solve
    claw.tfinal = 0.6
    status = claw.run()

    from pyclaw import plot
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
