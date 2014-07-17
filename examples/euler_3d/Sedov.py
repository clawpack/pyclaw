#!/usr/bin/env python
# encoding: utf-8

""" 
!-----------------------------------------------------------------------
! Description:  Test problem demonstrating a Sedov blast wave problem.
!   A step function "sphere" of energy is initialized at the center of
!   the domain.  This creates a shock that is evolved in time.
!
! Method:  This problem evolves the 3D Euler equations.
!   The primary variables are: 
!   density (rho), x,y, and z momentum (rho*u,rho*v,rho*w), and energy.
!-----------------------------------------------------------------------
"""
import numpy as np

# Constants
gamma = 1.4 # Ratio of Specific Heats
gamma1 = gamma - 1.
x0 = 0.0; y0 = 0.0; z0 = 0.0 # Sphere location
rmax = 0.10 # Radius of Sedov Sphere

# Setup Solver
def sedov(kernel_language='Fortran', solver_type='classic', use_petsc=False, splitting='unsplit', outdir='Sedov_output', output_format='hdf5', file_prefix='sedov', disable_output=False, mx=256, my=256, mz=256, tfinal=0.10, num_output_times=10):

    # Load PyClaw
    from clawpack import riemann
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    # Solver Settings
    if solver_type=='classic':
        solver = pyclaw.ClawSolver3D(riemann.euler_3D)
        solver.order = 2
        if splitting == "split":
            solver.dimensional_split = True
        else:
            solver.dimensional_split = False
            solver.transverse_waves = 22
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.num_ghost = 2
        solver.fwave = False
        solver.num_eqn = 5
        solver.num_waves = 3
        solver.cfl_max = 1.0
        solver.cfl_desired = 0.84
        solver.dt_initial = 3e-4
        solver.max_steps = 10000
    else:
        raise Exception('Unrecognized solver_type.')

    # Logging
    import logging
    solver.logger.setLevel(logging.DEBUG)

    # Initialize Domain
    x = pyclaw.Dimension('x', -1.0, 1.0, mx)
    y = pyclaw.Dimension('y', -1.0, 1.0, my)
    z = pyclaw.Dimension('z', -1.0, 1.0, mz)
    domain = pyclaw.Domain([x,y,z])

    # Define number of waves (eqn) and aux (eps,mu)
    num_aux = 0
    state = pyclaw.State(domain,solver.num_eqn, num_aux)
    state.problem_data['gamma']=gamma
    state.problem_data['gamma1']=gamma1
    
    # Initial Conditions
    grid = state.grid
    
    x = grid.x.centers
    y = grid.y.centers
    z = grid.z.centers
    X,Y,Z = np.meshgrid(x, y, z, indexing='ij')
    r = np.sqrt((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)

    state.q[0,:,:,:] = 1.0 # density (rho)
    state.q[1,:,:,:] = 0. # x-momentum (rho*u)
    state.q[2,:,:,:] = 0. # y-momentum (rho*v)
    state.q[3,:,:,:] = 0. # z-momentum (rho*w)
    
    p = np.zeros([mx,my,mz])
    vel_sqrd = np.zeros([mx,my,mz])
    eps = 1.0e-2
    Eblast = 0.851072
    ps = Eblast*gamma1/(4./3.*np.pi*rmax**3)
    p[:,:,:] = eps + ps*(r<=rmax) # pressure
    vel_sqrd[:,:,:] = (state.q[1,:,:,:]**2 \
                       + state.q[2,:,:,:]**2 \
                       + state.q[3,:,:,:]**2)/state.q[0,:,:,:]**2
    state.q[4,:,:,:] = p/gamma1 + 0.5*state.q[0,:,:,:]*vel_sqrd # energy (e)

    # Setup Boundary Conditions
    solver.all_bcs = pyclaw.BC.extrap

    # Solver Parameters
    claw = pyclaw.Controller()
    claw.verbosity = 4
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.output_format = output_format
    claw.output_file_prefix = file_prefix
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
