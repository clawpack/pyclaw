#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from scipy import integrate

gamma = 1.4
gamma1 = gamma - 1.
x0=0.5; y0=0.; r0=0.2
xshock = 0.2
pinf=5.

def inrad(y,x):
    return (np.sqrt((x-x0)**2+(y-y0)**2)<r0)

def ycirc(x,ymin,ymax):
    if r0**2>((x-x0)**2):
        return max(min(y0 + np.sqrt(r0**2-(x-x0)**2),ymax) - ymin,0.)
    else:
        return 0

def qinit(state,rhoin=0.1,bubble_shape='circle'):
    r"""
    Initialize data with a shock at x=xshock and a low-density bubble (of density rhoin)
    centersed at (x0,y0) with radius r0.
    """
    rhoout = 1.
    pout   = 1.
    pin    = 1.

    rinf = (gamma1 + pinf*(gamma+1.))/ ((gamma+1.) + gamma1*pinf)
    vinf = 1./np.sqrt(gamma) * (pinf - 1.) / np.sqrt(0.5*((gamma+1.)/gamma) * pinf+0.5*gamma1/gamma)
    einf = 0.5*rinf*vinf**2 + pinf/gamma1
    
    x =state.grid.x.centers
    y =state.grid.y.centers
    Y,X = np.meshgrid(y,x)
    if bubble_shape=='circle':
        r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    elif bubble_shape=='rectangle':
        z = np.dstack((np.abs(X-x0),np.abs(Y-y0)))
        r = np.max(z,axis=2)
    elif bubble_shape=='triangle':
        r = np.abs(X-x0) + np.abs(Y-y0)

    #First set the values for the cells that don't intersect the bubble boundary
    state.q[0,:,:] = rinf*(X<xshock) + rhoin*(r<=r0) + rhoout*(r>r0)*(X>xshock)
    state.q[1,:,:] = rinf*vinf*(X<xshock)
    state.q[2,:,:] = 0.
    state.q[3,:,:] = einf*(X<xshock) + (pin*(r<=r0) + pout*(r>r0)*(X>xshock))/gamma1
    state.q[4,:,:] = 1.*(r<=r0)

    #Now average for the cells on the edges of the bubble
    d2 = np.linalg.norm(state.grid.delta)/2.
    dx = state.grid.delta[0]
    dy = state.grid.delta[1]
    dx2 = state.grid.delta[0]/2.
    dy2 = state.grid.delta[1]/2.
    for i in xrange(state.q.shape[1]):
        for j in xrange(state.q.shape[2]):
            ydown = y[j]-dy2
            yup   = y[j]+dy2
            if abs(r[i,j]-r0)<d2:
                infrac,abserr = integrate.quad(ycirc,x[i]-dx2,x[i]+dx2,args=(ydown,yup),epsabs=1.e-8,epsrel=1.e-5)
                infrac=infrac/(dx*dy)
                state.q[0,i,j] = rhoin*infrac + rhoout*(1.-infrac)
                state.q[3,i,j] = (pin*infrac + pout*(1.-infrac))/gamma1
                state.q[4,i,j] = 1.*infrac


def auxinit(state):
    """
    aux[0,i,j] = y-coordinate of cell centers for cylindrical source terms
    """
    y=state.grid.y.centers
    for j,ycoord in enumerate(y):
        state.aux[0,:,j] = ycoord


def shockbc(state,dim,t,qbc,num_ghost):
    """
    Incoming shock at left boundary.
    """
    rinf = (gamma1 + pinf*(gamma+1.))/ ((gamma+1.) + gamma1*pinf)
    vinf = 1./np.sqrt(gamma) * (pinf - 1.) / np.sqrt(0.5*((gamma+1.)/gamma) * pinf+0.5*gamma1/gamma)
    einf = 0.5*rinf*vinf**2 + pinf/gamma1

    for i in xrange(num_ghost):
        qbc[0,i,...] = rinf
        qbc[1,i,...] = rinf*vinf
        qbc[2,i,...] = 0.
        qbc[3,i,...] = einf
        qbc[4,i,...] = 0.

def dq_Euler_radial(solver,state,dt):
    """
    Geometric source terms for Euler equations with radial symmetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a SharpClaw-style source term routine.
    """
    
    ndim = 2

    q   = state.q
    aux = state.aux

    rad = aux[0,:,:]

    rho = q[0,:,:]
    u   = q[1,:,:]/rho
    v   = q[2,:,:]/rho
    press  = gamma1 * (q[3,:,:] - 0.5*rho*(u**2 + v**2))

    dq = np.empty(q.shape)

    dq[0,:,:] = -dt*(ndim-1)/rad * q[2,:,:]
    dq[1,:,:] = -dt*(ndim-1)/rad * rho*u*v
    dq[2,:,:] = -dt*(ndim-1)/rad * rho*v*v
    dq[3,:,:] = -dt*(ndim-1)/rad * v * (q[3,:,:] + press)
    dq[4,:,:] = 0

    return dq

def step_Euler_radial(solver,state,dt):
    """
    Geometric source terms for Euler equations with radial symmetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine.
    """
    
    dt2 = dt/2.
    ndim = 2

    aux=state.aux
    q = state.q

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

    q[0,:,:] = q[0,:,:] - dt*(ndim-1)/rad * qstar[2,:,:]
    q[1,:,:] = q[1,:,:] - dt*(ndim-1)/rad * rho*u*v
    q[2,:,:] = q[2,:,:] - dt*(ndim-1)/rad * rho*v*v
    q[3,:,:] = q[3,:,:] - dt*(ndim-1)/rad * v * (qstar[3,:,:] + press)


def shockbubble(use_petsc=False,outdir='./_output',solver_type='classic'):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a bubble of dense gas that is impacted by a shock.
    """
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.euler_5wave_2D)
        solver.dq_src=dq_Euler_radial
        solver.weno_order=5
        solver.lim_type=2
    else:
        solver = pyclaw.ClawSolver2D(riemann.euler_5wave_2D)
        solver.dimensional_split = 0
        solver.transverse_waves = 2
        solver.limiters = [4,4,4,4,2]
        solver.step_source=step_Euler_radial

    solver.bc_lower[0]=pyclaw.BC.custom
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.wall
    solver.bc_upper[1]=pyclaw.BC.extrap

    #Aux variable in ghost cells doesn't matter
    solver.aux_bc_lower[0]=pyclaw.BC.extrap
    solver.aux_bc_upper[0]=pyclaw.BC.extrap
    solver.aux_bc_lower[1]=pyclaw.BC.extrap
    solver.aux_bc_upper[1]=pyclaw.BC.extrap

    # Initialize domain
    mx=160; my=40
    x = pyclaw.Dimension('x',0.0,2.0,mx)
    y = pyclaw.Dimension('y',0.0,0.5,my)
    domain = pyclaw.Domain([x,y])
    num_eqn = 5
    num_aux=1
    state = pyclaw.State(domain,num_eqn,num_aux)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1

    qinit(state)
    auxinit(state)

    solver.user_bc_lower=shockbc

    claw = pyclaw.Controller()
    claw.tfinal = 0.75
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.outdir = outdir

    return claw


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    from clawpack.pyclaw.examples import shock_bubble_interaction
    output = run_app_from_main(shockbubble,shock_bubble_interaction.setplot)
