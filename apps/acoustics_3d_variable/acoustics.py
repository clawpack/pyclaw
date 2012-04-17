#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def acoustics3D(iplot=False,htmlplot=False,use_petsc=False,outdir='./_output',solver_type='classic',**kwargs):
    """
    Example python script for solving the 3d acoustics equations.
    """

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver3D()
    else:
        raise Exception('Unrecognized solver_type.')

    from clawpack import riemann
    solver.rp = riemann.rp3_vc_acoustics
    solver.num_waves = 2
    solver.limiters = pyclaw.limiters.tvd.MC

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

    if 'test' in kwargs:
        test = kwargs['test']
        if test == 'homogeneous':
            app = 'test_homogeneous'
            solver.dimensional_split=True
            mx=256; my=4; mz=4
            zr = 1.0  # Impedance in right half
            cr = 1.0  # Sound speed in right half
        elif test == 'heterogeneous':
            app = 'test_heterogeneous'
            solver.dimensional_split=False
            solver.bc_lower[0]    =pyclaw.BC.wall
            solver.bc_lower[1]    =pyclaw.BC.wall
            solver.bc_lower[2]    =pyclaw.BC.wall
            solver.aux_bc_lower[0]=pyclaw.BC.wall
            solver.aux_bc_lower[1]=pyclaw.BC.wall
            solver.aux_bc_lower[2]=pyclaw.BC.wall
            mx=30; my=30; mz=30
            zr = 2.0  # Impedance in right half
            cr = 2.0  # Sound speed in right half

        else: raise Exception('Unrecognized test')
    else:
        app = 'example'
        solver.dimensional_split=False
        solver.limiters = pyclaw.limiters.tvd.MC

        solver.bc_lower[0]=pyclaw.BC.wall
        solver.bc_upper[0]=pyclaw.BC.periodic
        solver.bc_lower[1]=pyclaw.BC.wall
        solver.bc_upper[1]=pyclaw.BC.periodic
        solver.bc_lower[2]=pyclaw.BC.wall
        solver.bc_upper[2]=pyclaw.BC.periodic

        solver.aux_bc_lower[0]=pyclaw.BC.wall
        solver.aux_bc_upper[0]=pyclaw.BC.periodic
        solver.aux_bc_lower[1]=pyclaw.BC.wall
        solver.aux_bc_upper[1]=pyclaw.BC.periodic
        solver.aux_bc_lower[2]=pyclaw.BC.wall
        solver.aux_bc_upper[2]=pyclaw.BC.periodic
        mx=30; my=30; mz=30

    # Initialize domain
    x = pyclaw.Dimension('x',-1.0,1.0,mx)
    y = pyclaw.Dimension('y',-1.0,1.0,my)
    z = pyclaw.Dimension('z',-1.0,1.0,mz)
    domain = pyclaw.Domain([x,y,z])

    num_eqn = 4
    num_aux = 2 # density, sound speed
    state = pyclaw.State(domain,num_eqn,num_aux)

    zl = 1.0  # Impedance in left half
    cl = 1.0  # Sound speed in left half
    zr = 2.0  # Impedance in right half
    cr = 2.0  # Sound speed in right half

    grid = state.grid
    grid.compute_c_centers()
    X,Y,Z = grid._c_centers

    state.aux[0,:,:,:] = zl*(X<0.) + zr*(X>=0.) # Impedance
    state.aux[1,:,:,:] = cl*(X<0.) + cr*(X>=0.) # Sound speed

    x0 = -0.5; y0 = 0.; z0 = 0.
    if app == 'test_homogeneous':
        r = np.sqrt((X-x0)**2)
        width=0.2
        state.q[0,:,:,:] = (np.abs(r)<=width)*(1.+np.cos(np.pi*(r)/width))

    elif app == 'test_heterogeneous' or app == 'example':
        r = np.sqrt((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)
        width=0.1
        state.q[0,:,:,:] = (np.abs(r-0.3)<=width)*(1.+np.cos(np.pi*(r-0.3)/width))

    else: raise Exception('Unexpected application')
        
    state.q[1,:,:,:] = 0.
    state.q[2,:,:,:] = 0.
    state.q[3,:,:,:] = 0.

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir=outdir

    # Solve
    claw.tfinal = 2.0
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir,file_format=claw.output_format)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir,file_format=claw.output_format)


    pinitial=claw.frames[0].state.q[0,:,:,:].reshape(-1)
    pmiddle  =claw.frames[3].state.q[0,:,:,:].reshape(-1)
    pfinal  =claw.frames[claw.num_output_times].state.q[0,:,:,:].reshape(-1)

    #import matplotlib.pyplot as plt
    #for i in range(claw.num_output_times):
    #    plt.pcolor(claw.frames[i].state.q[0,:,:,mz/2])
    #    plt.figure()
    #plt.show()

    final_difference =np.prod(grid.delta)*np.linalg.norm(pfinal-pinitial,ord=1)
    middle_difference=np.prod(grid.delta)*np.linalg.norm(pmiddle-pinitial,ord=1)

    if app == 'example':
        print 'Final error: ', final_difference
        print 'Middle error: ', middle_difference

    return pfinal, final_difference

if __name__=="__main__":
    import sys
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics3D)

def test_3d_acoustics():
    """ tests against known classic results """


    def acoustics_verify_homogeneous(return_tuple):
        from clawpack.pyclaw.util import check_diff

        test_final_difference = return_tuple[1]
        return check_diff(0, test_final_difference, abstol=1e-4)

    def acoustics_verify_heterogeneous(return_tuple):
        import os
        import numpy as np
        from clawpack.pyclaw.util import check_diff

        test_pfinal = return_tuple[0]
        thisdir = os.path.dirname(__file__)
        verify_pfinal = np.loadtxt(os.path.join(thisdir,'verify_classic_heterogeneous.txt'))
        norm_err = np.linalg.norm(test_pfinal-verify_pfinal)
        return check_diff(0, norm_err, abstol=1e-4)

    from clawpack.pyclaw.util import gen_variants

    homogeneous_tests   = gen_variants(acoustics3D, acoustics_verify_homogeneous,
                                       python_kernel=False, solver_type='classic', test='homogeneous')

    heterogeneous_tests = gen_variants(acoustics3D, acoustics_verify_heterogeneous,
                                       python_kernel=False, solver_type='classic', test = 'heterogeneous')
    
    from itertools import chain
    for test in chain(homogeneous_tests, heterogeneous_tests):
        yield test
