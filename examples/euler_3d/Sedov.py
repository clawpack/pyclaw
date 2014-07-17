#!/usr/bin/env python
# encoding: utf-8

""" 
Test problem demonstrating a Sedov blast wave problem.
A spherical step function energy perturbation is initialized at the center of
the domain.  This creates an expanding shock wave.

This problem evolves the 3D Euler equations.
The primary variables are: 
    density (rho), x,y, and z momentum (rho*u,rho*v,rho*w), and energy.
"""
import numpy as np
from scipy import integrate

gamma = 1.4 # Ratio of Specific Heats
gamma1 = gamma - 1.
x0 = 0.0; y0 = 0.0; z0 = 0.0 # Sphere location
rmax = 0.10 # Radius of Sedov Sphere

def sphere_top(y, x):
    z2 = rmax**2 - (x-x0)**2 - (y-y0)**2
    if z2 < 0:
        return 0
    else:
        return np.sqrt(z2)

def sphere_bottom(y, x):
    return -sphere_top(y,x)

def f(y, x, zdown, zup):
    top = min(sphere_top(y,x), zup)
    bottom = min(top,max(sphere_bottom(y,x), zdown))
    return top-bottom


def sedov(kernel_language='Fortran', solver_type='classic', use_petsc=False,
          dimensional_split=False, outdir='Sedov_output', output_format='hdf5',
          disable_output=False, mx=64, my=64, mz=64,
          tfinal=0.10, num_output_times=10):

    from clawpack import riemann
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver3D(riemann.euler_3D)
        solver.dimensional_split = dimensional_split
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.cfl_max = 1.0
        solver.cfl_desired = 0.84
        solver.dt_initial = 3e-4
    else:
        raise Exception('Unrecognized solver_type.')

    x = pyclaw.Dimension('x', -1.0, 1.0, mx)
    y = pyclaw.Dimension('y', -1.0, 1.0, my)
    z = pyclaw.Dimension('z', -1.0, 1.0, mz)
    domain = pyclaw.Domain([x,y,z])

    state = pyclaw.State(domain,solver.num_eqn)

    state.problem_data['gamma']=gamma
    state.problem_data['gamma1']=gamma1
    
    grid = state.grid
    X,Y,Z = grid.p_centers
    r = np.sqrt((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)

    state.q[0,:,:,:] = 1.0 # density (rho)
    state.q[1,:,:,:] = 0. # x-momentum (rho*u)
    state.q[2,:,:,:] = 0. # y-momentum (rho*v)
    state.q[3,:,:,:] = 0. # z-momentum (rho*w)
    
    background_pressure = 1.0e-2
    Eblast = 0.851072
    pressure_in = Eblast*gamma1/(4./3.*np.pi*rmax**3)
    state.q[4,:,:,:] = background_pressure/gamma1 # energy (e)

    # Compute cell fraction inside initial perturbed sphere
    dx, dy, dz = state.grid.delta
    dx2, dy2, dz2 = [d/2. for d in state.grid.delta]
    dmax = max(state.grid.delta)

    for i in xrange(state.q.shape[1]):
        for j in xrange(state.q.shape[2]):
            for k in xrange(state.q.shape[3]):
                if r[i,j,k] - dmax > rmax:
                    continue
                xdown = X[i,j,k] - dx2
                xup   = X[i,j,k] + dx2
                ydown = lambda x : Y[i,j,k] - dy2
                yup   = lambda x : Y[i,j,k] + dy2
                zdown = Z[i,j,k] - dz2
                zup   = Z[i,j,k] + dz2

                infrac,abserr = integrate.dblquad(f,xdown,xup,ydown,yup,args=(zdown,zup),epsabs=1.e-3,epsrel=1.e-2)
                infrac=infrac/(dx*dy*dz)

                p = background_pressure + pressure_in*infrac # pressure
                state.q[4,i,j,k] = p/gamma1 # energy (e)

    solver.all_bcs = pyclaw.BC.extrap

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
    output = run_app_from_main(sedov)
