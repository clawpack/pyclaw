#!/usr/bin/env python
# encoding: utf-8
r"""
Three-dimensional variable-coefficient acoustics
==============================================

Solve the variable-coefficient acoustics equations in 3D:

.. math:: 
    p_t + K(x,y,z) (u_x + v_y + w_z) & = 0 \\ 
    u_t + p_x / \rho(x,y,z) & = 0 \\
    v_t + p_y / \rho(x,y,z) & = 0 \\
    w_t + p_z / \rho(x,y,z) & = 0 \\

Here p is the pressure, (u,v,w) is the velocity, :math:`K(x,y,z)` is the bulk modulus,
and :math:`\rho(x,y,z)` is the density.

This example shows how to solve a problem with variable coefficients.
The left and right halves of the domain consist of different materials.
"""
 
import numpy as np

def setup(use_petsc=False,outdir='./_output',solver_type='classic',
          mx=30,my=30,mz=30,disable_output=False,problem='heterogeneous',**kwargs):
    """
    Example python script for solving the 3d acoustics equations.
    """
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver3D(riemann.vc_acoustics_3D)
        solver.limiters = pyclaw.limiters.tvd.MC
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver3D(riemann.vc_acoustics_3D)
        
    else:
        raise Exception('Unrecognized solver_type.')


    solver.bc_lower[0]=pyclaw.BC.periodic
    solver.bc_upper[0]=pyclaw.BC.periodic
    solver.bc_lower[1]=pyclaw.BC.periodic
    solver.bc_upper[1]=pyclaw.BC.periodic
    solver.bc_lower[2]=pyclaw.BC.periodic
    solver.bc_upper[2]=pyclaw.BC.periodic

    solver.aux_bc_lower[0]=pyclaw.BC.periodic
    solver.aux_bc_upper[0]=pyclaw.BC.periodic
    solver.aux_bc_lower[1]=pyclaw.BC.periodic
    solver.aux_bc_upper[1]=pyclaw.BC.periodic
    solver.aux_bc_lower[2]=pyclaw.BC.periodic
    solver.aux_bc_upper[2]=pyclaw.BC.periodic

    zl = 1.0  # Impedance in left half
    cl = 1.0  # Sound speed in left half

    if problem == 'homogeneous':
        if solver_type=='classic':
            solver.dimensional_split=True
        else:
            solver.lim_type = 1

        solver.limiters = [4]
        
        mx=mx; my=my; mz=mz # Grid resolution

        zr = 1.0  # Impedance in right half
        cr = 1.0  # Sound speed in right half

    if problem == 'heterogeneous':
        if solver_type=='classic':
            solver.dimensional_split=False
        
        solver.bc_lower[0]    =pyclaw.BC.wall
        solver.bc_lower[1]    =pyclaw.BC.wall
        solver.bc_lower[2]    =pyclaw.BC.wall
        solver.aux_bc_lower[0]=pyclaw.BC.wall
        solver.aux_bc_lower[1]=pyclaw.BC.wall
        solver.aux_bc_lower[2]=pyclaw.BC.wall

        mx=mx; my=my; mz=mz # Grid resolution

        zr = 2.0  # Impedance in right half
        cr = 2.0  # Sound speed in right half

    solver.limiters = pyclaw.limiters.tvd.MC

    # Initialize domain
    x = pyclaw.Dimension(-1.0,1.0,mx,name='x')
    y = pyclaw.Dimension(-1.0,1.0,my,name='y')
    z = pyclaw.Dimension(-1.0,1.0,mz,name='z')
    domain = pyclaw.Domain([x,y,z])

    num_eqn = 4
    num_aux = 2 # density, sound speed
    state = pyclaw.State(domain,num_eqn,num_aux)

    X,Y,Z = state.grid.p_centers

    state.aux[0,:,:,:] = zl*(X<0.) + zr*(X>=0.) # Impedance
    state.aux[1,:,:,:] = cl*(X<0.) + cr*(X>=0.) # Sound speed

    # Set initial density
    x0 = -0.5; y0 = 0.; z0 = 0.
    if problem == 'homogeneous':
        r = np.sqrt((X-x0)**2)
        width=0.2
        state.q[0,:,:,:] = (np.abs(r)<=width)*(1.+np.cos(np.pi*(r)/width))
    elif problem == 'heterogeneous':
        r = np.sqrt((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)
        width=0.1
        state.q[0,:,:,:] = (np.abs(r-0.3)<=width)*(1.+np.cos(np.pi*(r-0.3)/width))
    else: 
        raise Exception('Unrecognized problem name')
        
    # Set initial velocities to zero
    state.q[1,:,:,:] = 0.
    state.q[2,:,:,:] = 0.
    state.q[3,:,:,:] = 0.

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
       claw.output_format = None
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = 2.0

    return claw


if __name__=="__main__":
    import sys
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
