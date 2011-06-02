#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from petsc4py import PETSc



def qinit(grid):

    # Set initial conditions for q.
    # Sample scalar equation with data that is piecewise constant with
    # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    #     0.1  otherwise

    grid.zeros_q()
    
    x =grid.x.center
    y =grid.y.center
    #q2=np.arange(len(q.reshape([-1])))
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] > 0.0 and x[i] < 0.5 and y[j]>0.0 and y[j] < 0.5:
                grid.q[:,i,j] = 1.0
            else:
                grid.q[:,i,j] = 0.1
                
def advection2D(iplot=False,use_petsc=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    """
    Example python script for solving the 2d advection equation.
    """

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()

    solver.mthbc_lower[0] = pyclaw.BC.periodic
    solver.mthbc_upper[0] = pyclaw.BC.periodic
    solver.mthbc_lower[1] = pyclaw.BC.periodic
    solver.mthbc_upper[1] = pyclaw.BC.periodic

    mx=80; my=80
    # Initialize grids and solutions
    x = pyclaw.Dimension('x',0.0,1.0,mx)
    y = pyclaw.Dimension('y',0.0,1.0,my)
    grid = pyclaw.Grid([x,y])

    grid.mbc = solver.mbc

    grid.aux_global['u']=0.6
    grid.aux_global['v']=0.4

    grid.meqn = 1
    qinit(grid)

    solver.dim_split = 1
    solver.cfl_max=0.5
    solver.cfl_desired = 0.45
    solver.mwaves=1
    solver.mthlim = pyclaw.limiters.vanleer

    claw = pyclaw.Controller()
    claw.tfinal = 1.0
    claw.solution = pyclaw.Solution(grid)
    claw.solver = solver
    claw.outdir = outdir

    #Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.plotHTML(outdir=outdir)
    if iplot:     pyclaw.plot.plotInteractive(outdir=outdir)


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        advection2D(*args,**kwargs)
    else: advection2D()
