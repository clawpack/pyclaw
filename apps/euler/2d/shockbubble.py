#!/usr/bin/env python
# encoding: utf-8

import numpy as np

gamma = 1.4
gamma1 = gamma - 1.

def qinit(state,x0=0.5,y0=0.,r0=0.2,rhoin=0.1,pinf=5.):
    rhoout = 1.
    pout   = 1.
    pin    = 1.

    rinf = (gamma1 + pinf*(gamma+1.))/ ((gamma+1.) + gamma1*pinf)
    vinf = 1./np.sqrt(gamma) * (pinf - 1.) / np.sqrt(0.5*((gamma+1.)/gamma) * pinf+0.5*gamma1/gamma)
    einf = 0.5*rinf*vinf**2 + pinf/gamma1
    
    x =state.grid.x.center
    y =state.grid.y.center
    Y,X = np.meshgrid(y,x)
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    
    state.q[0,:,:] = rhoin*(r<=r0) + rhoout*(r>r0)
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.
    state.q[3,:,:] = (pin*(r<=r0) + pout*(r>r0))/gamma1
    state.q[4,:,:] = 1.*(r<=r0)


def auxinit(state):
    """
    aux[0,i,j] = y-coordinate of cell center for cylindrical source terms
    """
    state.maux=1
    x=state.grid.x.center
    y=state.grid.y.center
    for j,ycoord in enumerate(y):
        state.aux[0,:,j] = ycoord


def shockbc(grid,dim,t,qbc,mbc):
    """
    Incoming shock at left boundary.
    """
    if dim.nstart == 0:

        pinf=5.
        rinf = (gamma1 + pinf*(gamma+1.))/ ((gamma+1.) + gamma1*pinf)
        vinf = 1./np.sqrt(gamma) * (pinf - 1.) / np.sqrt(0.5*((gamma+1.)/gamma) * pinf+0.5*gamma1/gamma)
        einf = 0.5*rinf*vinf**2 + pinf/gamma1

        for i in xrange(mbc):
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
    srcterm_language = 'Python'
    state = solutions['n'].states[0]
    aux=state.aux
    grid = state.grid
    q = state.q
    
    if  srcterm_language == 'Python':
        dt2 = dt/2.
        press = 0.
        ndim = 2


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
    else:
        from classic2 import src2
        qbc = solver.append_ghost_cells(state)
        auxbc = solver.append_ghost_cells_to_aux(state)
        maxmx, maxmy, mx, my = grid.ng[0], grid.ng[1], grid.ng[0], grid.ng[1]
        mbc = solver.mbc
        xlower, ylower, dx, dy =0,0,0,0 # it is not used
        #t, dt = passed
        qbc = src2(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,qbc,auxbc,t,dt)
        solver.set_global_q(state, qbc)

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

    #Aux variable in ghost cells doesn't matter
    solver.mthauxbc_lower[0]=pyclaw.BC.outflow
    solver.mthauxbc_upper[0]=pyclaw.BC.outflow
    solver.mthauxbc_lower[1]=pyclaw.BC.outflow
    solver.mthauxbc_upper[1]=pyclaw.BC.outflow

    # Initialize grid
    mx=600; my=200
    x = pyclaw.Dimension('x',0.0,2.0,mx)
    y = pyclaw.Dimension('y',0.0,0.5,my)
    grid = pyclaw.Grid([x,y])
    state = pyclaw.State(grid)

    state.aux_global['gamma']= gamma
    state.aux_global['gamma1']= gamma1
    state.meqn = 5

    qinit(state)
    auxinit(state)

    solver.dim_split = 0
    solver.limiters = [4,4,4,4,2]
    solver.user_bc_lower=shockbc
    solver.src=euler_rad_src
    solver.src_split = 1

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = 0.4
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.nout = 10
    claw.outdir = outdir

    # Solve
    status = claw.run()

    if htmlplot:  pyclaw.plot.plotHTML(outdir=outdir)
    if iplot:     pyclaw.plot.plotInteractive(outdir=outdir)

    return claw.solution.q

if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        shockbubble(*args,**kwargs)
    else: shockbubble()
