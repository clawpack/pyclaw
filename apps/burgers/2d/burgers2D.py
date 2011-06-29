#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from petsc4py import PETSc
from petclaw.grid import Dimension
from petclaw.grid import Grid
from pyclaw.solution import Solution
from petclaw.clawpack import PetClawSolver2D
from pyclaw.controller import Controller
import math
from petclaw import plot


def qinit(grid):
    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)

    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    q[0,:,:] = 0.1 + 0.9* (X>0.1)*(X<0.6)*(Y>0.1)*(Y<0.6)
    grid.q=q


def burgers2D(iplot=False,petscPlot=True,useController=True,htmlplot=False,outdir='./_output'):
    """
    Example python script for solving the 2d Burgers' equations.
    """
    
    # Initialize grid
    mx=50; my=50
    x = Dimension('x',0.0,1.0,mx,mthbc_lower=2,mthbc_upper=2)
    y = Dimension('y',-0.0,1.0,my,mthbc_lower=2,mthbc_upper=2)
    grid = Grid([x,y])

    grid.meqn = 1
    grid.mbc = 2
    tfinal = 1.0
    qinit(grid)
    inital_solution = Solution(grid)

    solver = PetClawSolver2D()
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.4
    solver.mwaves = 1
    solver.mthlim = [4]

    claw = Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.output_format = 'petsc'
    claw.tfinal = tfinal
    claw.solutions['n'] = inital_solution
    claw.solver = solver
    claw.outdir = outdir

    # Solve
    status = claw.run()

    if htmlplot:  plot.html_plot(outdir=outdir,format=output_format)
    if iplot:     plot.interactive_plot(outdir=outdir,format=output_format)


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        burgers2D(*args,**kwargs)
    else: burgers2D()
