#!/usr/bin/env python
# encoding: utf-8
r"""
Variable-coefficient elasticity example.
"""
from __future__ import absolute_import
import numpy as np
from six.moves import range
t0wall  = 0.025
tperiod = 0.05

def moving_wall_bc(state,dim,t,qbc,num_ghost):
    x = state.grid.x.centers_with_ghost(num_ghost)[:num_ghost]
    if t<t0wall:
        s = np.sin(np.pi*t/tperiod)
    else:
        s = 0.

    for i in range(num_ghost):
        # First reflect-extrapolate
        qbc[:,i,:] = qbc[:,2*num_ghost-i-1,:]
        # Now set velocity
        qbc[3,i,:] = 2.0*s - qbc[3,i,:]
        qbc[4,i,:] =       - qbc[4,i,:]
    

#def no_stress_bc(state,dim,t,qbc,num_ghost):
#    """No-stress boundary condition: sigma_{12} = sigma_{11} = 0"""
#    if state.grid.on_lower_boundary[idim]:
#        jghost = 
#    # First extrapolate
#    tmp = qbc[:,:,self.num_ghost]
#    tmp = np.tile(tmp,(1,1,num_ghost))
#    qbc[:,i,: = tmp
#
#    # Then negate the sig12 and sig11 components

def integrate_displacement(solver,state):
    aux[5,:,:] = aux[5,:,:] + dt*q[3,:,:]
    aux[6,:,:] = aux[6,:,:] + dt*q[4,:,:]


def inclusion():
    from clawpack import pyclaw
    from clawpack import riemann

    solver=pyclaw.ClawSolver2D(riemann.vc_elasticity_2D)
    solver.dimensional_split = False
    solver.transverse_waves = 2
    solver.limiters = pyclaw.limiters.tvd.MC


    mx = 200
    my = 100
    num_aux = 7
    domain = pyclaw.Domain( (0.,0.),(2.,1.),(mx,my) )
    state = pyclaw.State(domain,solver.num_eqn,num_aux)
    solution = pyclaw.Solution(state,domain)


    solver.bc_lower[0] = pyclaw.BC.custom
    solver.user_bc_lower=moving_wall_bc
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.periodic  # No stress
    solver.bc_upper[1] = pyclaw.BC.periodic  # No stress

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.extrap
    solver.aux_bc_upper[1] = pyclaw.BC.extrap

    rho1 = 1.0
    lam1 = 200.
    mu1  = 100.

    rho2 = 1.0
    lam2 = 2.0
    mu2  = 1.0


    # set aux arrays
    #  aux[0,i,j] = density rho in (i,j) cell
    #  aux[1,i,j] = lambda in (i,j) cell
    #  aux[2,i,j] = mu in (i,j) cell
    #  aux[3,i,j] = cp in (i,j) cell
    #  aux[4,i,j] = cs in (i,j) cell
    #  aux[5,i,j] = xdisp in (i,j) cell
    #  aux[6,i,j] = ydisp in (i,j) cell

    xx,yy = domain.grid.p_centers
    inbar = (0.5<xx)*(xx<1.5)*(0.4<yy)*(yy<0.6)
    outbar = 1 - inbar
    aux = state.aux
    aux[0,:,:] = rho1 * inbar + rho2 * outbar
    aux[1,:,:] = lam1 * inbar + lam2 * outbar
    aux[2,:,:] = mu1  * inbar + mu2  * outbar
    bulk       = aux[1,:,:] + 2.*aux[2,:,:]
    aux[3,:,:] = np.sqrt(bulk/aux[0,:,:])
    aux[4,:,:] = np.sqrt(aux[2,:,:]/aux[0,:,:])
    aux[5,:,:] = 0.
    aux[6,:,:] = 0.


    # set initial condition
    state.q[:,:,:] = 0.


    claw = pyclaw.Controller()
    claw.solver = solver
    claw.solution = solution
    claw.num_output_times = 20
    claw.tfinal = 0.5

    return claw


if __name__ == '__main__':
    claw = inclusion()
    claw.run()

    from clawpack.pyclaw import plot
    plot.interactive_plot()

