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

    grid.zeros_q()
    grid.q[0,:,:] = rhoin*(r<=r0) + rhoout*(r>r0)
    grid.q[1,:,:] = 0.
    grid.q[2,:,:] = 0.
    grid.q[3,:,:] = (pin*(r<=r0) + pout*(r>r0))/gamma1
    grid.q[4,:,:] = 1.*(r<=r0)

def auxinit(grid):
    """
    aux[0,i,j] = y-coordinate of cell center for cylindrical source terms
    """
    x=grid.x.centerghost
    y=grid.y.centerghost
    grid.zeros_aux(1)
    #aux=np.empty([1,len(x),len(y)], order='F')
    for j,ycoord in enumerate(y):
        grid.aux[0,:,j] = ycoord

def shockbc(grid,dim,t,qbc):
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


def shockbubble(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a bubble of dense gas that is impacted by a shock.
    """

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()
    else:
        solver = pyclaw.ClawSolver2D()

    solver.mwaves = 5
    solver.mthbc_lower[0]=pyclaw.BC.custom
    solver.mthbc_upper[0]=pyclaw.BC.outflow
    solver.mthbc_lower[1]=pyclaw.BC.reflecting
    solver.mthbc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=640; my=160
    x = pyclaw.Dimension('x',0.0,2.0,mx)
    y = pyclaw.Dimension('y',0.0,0.5,my)
    grid = pyclaw.Grid([x,y])

    grid.aux_global['gamma']= gamma
    grid.aux_global['gamma1']= gamma1
    grid.meqn = 5
    grid.mbc = solver.mbc

    qinit(grid)
    auxinit(grid)

    solver.dim_split = 0
    solver.mthlim = [4,4,4,4,2]
    solver.user_bc_lower=shockbc
    solver.src=euler_rad_src
    solver.src_split = 1

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = 0.75
    claw.solution = pyclaw.Solution(grid)
    claw.solver = solver
    claw.nout = 10
    claw.outdir = outdir

    # Solve
    status = claw.run()

    from pyclaw import plot
    if htmlplot:  plot.plotHTML(outdir=outdir)
    if iplot:     plot.plotInteractive(outdir=outdir)

    return claw.solution.q

if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        shockbubble(*args,**kwargs)
    else: shockbubble()
