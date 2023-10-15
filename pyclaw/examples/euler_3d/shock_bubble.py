#!/usr/bin/env python
# encoding: utf-8
r""" 
3D shock-bubble interaction problem.
A planar shock wave impacts a spherical region of low density.

This problem involves the 3D Euler equations:
.. math::
    \rho_t + (\rho u)_x + (\rho v)_y + (\rho w)_z & = 0 \\
    (\rho u)_t + (\rho u^2 + p)_x + (\rho uv)_y & = 0 \\
    (\rho v)_t + (\rho uv)_x + (\rho v^2 + p)_y + (\rho vw)_z & = 0 \\
    (\rho w)_t + (\rho uw)_x + (\rho vw)_y + (\rho w^2 + p)_z & = 0 \\
    E_t + \nabla \cdot (u (E + p) ) & = 0.


The conserved quantities are: 
    density (rho), x-,y-, and z-momentum (rho*u,rho*v,rho*w), and energy.
"""
from __future__ import absolute_import
import numpy as np
from scipy import integrate
from six.moves import range

gamma = 1.4 # Ratio of Specific Heats
gamma1 = gamma - 1.
x0 = 0.5; y0 = 0.; z0 = 0. # Bubble location
r_bubble = 0.2 # Bubble radius

# Ambient state
rhoout = 1.0
pout   = 1.0

# Bubble state
rhoin  = 0.1
pin    = 1.0

xshock = 0.2 # Initial shock wave location

# State behind shock wave
p_shock = 5.0
rho_shock = (gamma1 + p_shock*(gamma+1.))/ ((gamma+1.) + gamma1*p_shock)
v_shock = (p_shock - 1.) / np.sqrt(0.5 * ((gamma+1.) * p_shock+gamma1))
e_shock = 0.5*rho_shock*v_shock**2 + p_shock/gamma1


def bubble(y, x, zdown, zup, which):
    "Used to compute how much of each cell is in the bubble."
    def sphere_top(y, x, which):
        z2 = r_bubble**2 - (x-x0)**2 - (y-y0)**2
        if z2 < 0:
            return 0
        else:
            return z0 + np.sqrt(z2)

    def sphere_bottom(y, x, which):
        z2 = (r_bubble**2 - (x-x0)**2 - (y-y0)**2)
        if z2 < 0:
            return 0
        else:
            return z0 - np.sqrt(z2)

    top = min(sphere_top(y,x, which), zup)
    bottom = min(top,max(sphere_bottom(y,x, which), zdown))
    return top-bottom

def incoming_shock(state,dim,t,qbc,auxbc,num_ghost):
    """
    Incoming shock at x=0 boundary.
    """
    for i in range(num_ghost):
        qbc[0,i,...] = rho_shock
        qbc[1,i,...] = rho_shock*v_shock
        qbc[2,i,...] = 0.
        qbc[3,i,...] = 0.
        qbc[4,i,...] = e_shock


def setup(kernel_language='Fortran', solver_type='classic', use_petsc=False,
          dimensional_split=False, outdir='_output', output_format='hdf5',
          disable_output=False, num_cells=(256,64,64),
          tfinal=0.6, num_output_times=10):

    from clawpack import riemann
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver3D(riemann.euler_3D)
        solver.dimensional_split = dimensional_split
    else:
        raise Exception('Unrecognized solver_type.')

    x = pyclaw.Dimension(0.0, 2.0, num_cells[0], name='x')
    y = pyclaw.Dimension(0.0, 0.5, num_cells[1], name='y')
    z = pyclaw.Dimension(0.0, 0.5, num_cells[2], name='z')
    domain = pyclaw.Domain([x,y,z])

    solver.all_bcs = pyclaw.BC.extrap
    solver.bc_lower[0]  = pyclaw.BC.custom
    solver.user_bc_lower = incoming_shock
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_lower[2] = pyclaw.BC.wall

    state = pyclaw.State(domain,solver.num_eqn)

    state.problem_data['gamma']  = gamma
    
    grid = state.grid
    X,Y,Z = grid.p_centers
    r0 = np.sqrt((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)

    state.q[0,:,:,:] = rho_shock*(X<xshock) + rhoout*(X>=xshock) # density (rho)
    state.q[1,:,:,:] = rho_shock*v_shock*(X<xshock) # x-momentum (rho*u)
    state.q[2,:,:,:] = 0. # y-momentum (rho*v)
    state.q[3,:,:,:] = 0. # z-momentum (rho*w)
    state.q[4,:,:,:] = e_shock*(X<xshock) + pout/gamma1*(X>=xshock) # energy (e)

    # Compute cell fraction inside bubble
    dx, dy, dz = state.grid.delta
    dx2, dy2, dz2 = [d/2. for d in state.grid.delta]
    dmax = max(state.grid.delta)

    for i in range(state.q.shape[1]):
        for j in range(state.q.shape[2]):
            for k in range(state.q.shape[3]):
                if (r0[i,j,k] - dmax > r_bubble):
                    continue
                xdown = X[i,j,k] - dx2
                xup   = X[i,j,k] + dx2
                ydown = lambda x : Y[i,j,k] - dy2
                yup   = lambda x : Y[i,j,k] + dy2
                zdown = Z[i,j,k] - dz2
                zup   = Z[i,j,k] + dz2

                infrac,abserr = integrate.dblquad(bubble,xdown,xup,ydown,yup,
                                                  args=(zdown,zup,0),
                                                  epsabs=1.e-3,epsrel=1.e-2)
                infrac=infrac/(dx*dy*dz)

                state.q[0,i,j,k] = rhoout*(1.-infrac) + rhoin*infrac
                state.q[4,i,j,k] = (pout*(1.-infrac) + pin*infrac)/gamma1 # energy (e)

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.output_format = output_format
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.tfinal = tfinal
    claw.num_output_times = num_output_times
    claw.outdir = outdir

    return claw

# __main__()
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
