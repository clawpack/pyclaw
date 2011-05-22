#!/usr/bin/env python
# encoding: utf-8

import numpy as np

gamma = 1.4
gamma1 = gamma - 1.

def qinit(grid,x0=0.5,y0=0.,r0=0.2,rhoin=0.1,pinf=5.):
    rhoout = 1.
    pout   = 1.
    pin    = 1.

    rinf = (gamma1 + pinf*(gamma+1.))/ ((gamma+1.) + gamma1*pinf)
    vinf = 1./np.sqrt(gamma) * (pinf - 1.) / np.sqrt(0.5*((gamma+1.)/gamma) * pinf+0.5*gamma1/gamma)
    einf = 0.5*rinf*vinf**2 + pinf/gamma1
    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)

    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    q[0,:,:] = rhoin*(r<=r0) + rhoout*(r>r0)
    q[1,:,:] = 0.
    q[2,:,:] = 0.
    q[3,:,:] = (pin*(r<=r0) + pout*(r>r0))/gamma1
    q[4,:,:] = 1.*(r<=r0)
    grid.q=q

def auxinit(grid):
    """
    aux[0,i,j] = y-coordinate of cell center for cylindrical source terms
    """
    x=grid.x.centerghost
    y=grid.y.centerghost
    aux=np.empty([1,len(x),len(y)], order='F')
    for j,ycoord in enumerate(y):
        aux[0,:,j] = ycoord
    grid.aux=aux

def shockbc(grid,dim,qbc):
    """
    Incoming shock at left boundary.
    """
    if dim.nstart == 0:

        pinf=5.
        rinf = (gamma1 + pinf*(gamma+1.))/ ((gamma+1.) + gamma1*pinf)
        vinf = 1./np.sqrt(gamma) * (pinf - 1.) / np.sqrt(0.5*((gamma+1.)/gamma) * pinf+0.5*gamma1/gamma)
        einf = 0.5*rinf*vinf**2 + pinf/gamma1

        for i in xrange(grid.mbc):
            qbc[0,i,...] = rinf
            qbc[1,i,...] = rinf*vinf
            qbc[2,i,...] = 0.
            qbc[3,i,...] = einf
            qbc[4,i,...] = 0.

def euler_rad_src(solver,solutions,t,dt):
    """
    Geometric source terms for Euler equations with radial symmetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    """
    
    dt2 = dt/2.
    press = 0.
    ndim = 2
    mbc=solver.mbc

    aux=solutions['n'].grids[0].aux[:,mbc:-mbc,mbc:-mbc]
    q = solutions['n'].grids[0].q

    rad = aux[0,:,:]

    rho = q[0,:,:]
    u   = q[1,:,:]/rho
    v   = q[2,:,:]/rho
    press  = gamma1 * (q[3,:,:] - 0.5*rho*(u**2 + v**2))

    qstar = np.empty(q.shape)

    qstar[0,:,:] = q[0,:,:] - dt2*(ndim-1)/rad * q[2,:,:]
    qstar[1,:,:] = q[1,:,:] - dt2*(ndim-1)/rad * rho*u*v
    qstar[2,:,:] = q[2,:,:] - dt2*(ndim-1)/rad * rho*v*v
    qstar[3,:,:] = q[3,:,:] - dt2*(ndim-1)/rad * v * (q[3,:,:] + press)

    rho = qstar[0,:,:]
    u   = qstar[1,:,:]/rho
    v   = qstar[2,:,:]/rho
    press  = gamma1 * (qstar[3,:,:] - 0.5*rho*(u**2 + v**2))

    q[0,:,:] = q[0,:,:] - dt2*(ndim-1)/rad * qstar[2,:,:]
    q[1,:,:] = q[1,:,:] - dt2*(ndim-1)/rad * rho*u*v
    q[2,:,:] = q[2,:,:] - dt2*(ndim-1)/rad * rho*v*v
    q[3,:,:] = q[3,:,:] - dt2*(ndim-1)/rad * v * (qstar[3,:,:] + press)


def shockbubble(use_PETSc=False,iplot=False,htmlplot=False,outdir='./_output'):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a bubble of dense gas that is impacted by a shock.
    """

    if use_PETSc:
        from petclaw.grid import Dimension, Grid
        from petclaw.evolve.clawpack import PetClawSolver2D as mySolver
        output_format = 'petsc'
    else:
        from pyclaw.grid import Dimension, Grid
        from pyclaw.evolve.clawpack import ClawSolver2D as mySolver
        output_format = 'ascii'

    from pyclaw.solution import Solution
    from pyclaw.controller import Controller
    from petclaw import plot

    # Initialize grid
    mx=320; my=80
    x = Dimension('x',0.0,2.0,mx,mthbc_lower=0,mthbc_upper=1)
    y = Dimension('y',0.0,0.5,my,mthbc_lower=3,mthbc_upper=1)
    grid = Grid([x,y])

    grid.aux_global['gamma']= gamma
    grid.aux_global['gamma1']= gamma1
    from classic2 import cparam
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    grid.meqn = 5
    grid.mbc = 2
    tfinal = 0.75

    if use_PETSc:
        # Initialize petsc Structures for q
        grid.init_q_petsc_structures()

    qinit(grid)
    auxinit(grid)
    initial_solution = Solution(grid)

    solver = mySolver()
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 5
    solver.mthlim = [4,4,4,4,2]
    solver.dt=0.005
    solver.user_bc_lower=shockbc
    solver.src=euler_rad_src
    solver.src_split = 1

    claw = Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.output_format = output_format
    claw.tfinal = tfinal
    claw.solutions['n'] = initial_solution
    claw.solver = solver
    claw.nout = 10
    claw.outdir = outdir

    # Solve
    status = claw.run()

    if htmlplot:  plot.plotHTML(outdir=outdir,format=output_format)
    if iplot:     plot.plotInteractive(outdir=outdir,format=output_format)


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        shockbubble(*args,**kwargs)
    else: shockbubble()
