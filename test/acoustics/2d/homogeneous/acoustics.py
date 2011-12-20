#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(state,width=0.2):
    
    grid = state.grid
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt(X**2 + Y**2)

    state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.


def acoustics2D(use_petsc=False,kernel_language='Fortran',iplot=False,htmlplot=False,solver_type='classic', outdir = './_output', nout = 10):
    """
    Example python script for solving the 2d acoustics equations.
    """
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D()
        solver.dim_split = 1
        solver.mwaves = 2
        solver.limiters = [4]*solver.mwaves
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()
        solver.mwaves = 2

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.bc_lower[0] = pyclaw.BC.outflow
    solver.bc_upper[0] = pyclaw.BC.outflow
    solver.bc_lower[1] = pyclaw.BC.outflow
    solver.bc_upper[1] = pyclaw.BC.outflow

    # Initialize grid
    mx=100; my=100
    x = pyclaw.grid.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.grid.Dimension('y',-1.0,1.0,my)
    grid = pyclaw.grid.Grid([x,y])
    meqn = 3
    state = pyclaw.State(grid,meqn)

    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    state.aux_global['rho']= rho
    state.aux_global['bulk']=bulk
    state.aux_global['zz']= zz
    state.aux_global['cc']=cc

    tfinal = 0.12

    qinit(state)
    initial_solution = pyclaw.Solution(state)

    solver.dt_initial=np.min(grid.d)/state.aux_global['cc']*solver.cfl_desired

    claw = pyclaw.Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.tfinal = tfinal
    claw.solution = initial_solution
    claw.solver = solver
    claw.outdir = outdir
    claw.nout = nout

    # Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot()
    if iplot:     pyclaw.plot.interactive_plot()

    if use_petsc:
        pressure=claw.frames[claw.nout].state.gqVec.getArray().reshape([state.meqn,grid.ng[0],grid.ng[1]],order='F')[0,:,:]
    else:
        pressure=claw.frames[claw.nout].state.q[0,:,:]
    return pressure


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        acoustics2D(*args,**kwargs)
    else: acoustics2D()
