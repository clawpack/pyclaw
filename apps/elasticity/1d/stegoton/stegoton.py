#!/usr/bin/env python
# encoding: utf-8
import numpy as np

def qinit(state,ic=2,a2=1.0,xupper=600.):
    x = state.grid.x.center
    
    if ic==1: #Zero ic
        state.q[:,:] = 0.
    elif ic==2:
        # Gaussian
        sigma = a2*np.exp(-((x-xupper/2.)/10.)**2.)
        state.q[0,:] = np.log(sigma+1.)/state.aux[1,:]
        state.q[1,:] = 0.


def setaux(x,rhoB=4,KB=4,rhoA=1,KA=1,alpha=0.5,xlower=0.,xupper=600.,bc=2):
    aux = np.empty([3,len(x)],order='F')
    xfrac = x-np.floor(x)
    #Density:
    aux[0,:] = rhoA*(xfrac<alpha)+rhoB*(xfrac>=alpha)
    #Bulk modulus:
    aux[1,:] = KA  *(xfrac<alpha)+KB  *(xfrac>=alpha)
    aux[2,:] = 0. # not used
    return aux

    
def b4step(solver,solution):
    #Reverse velocity at trtime
    #Note that trtime should be an output point
    state = solution.states[0]
    if state.t>=state.aux_global['trtime']-1.e-10 and not state.aux_global['trdone']:
        #print 'Time reversing'
        state.q[1,:]=-state.q[1,:]
        solution.state.q=state.q
        state.aux_global['trdone']=True
        if state.t>state.aux_global['trtime']:
            print 'WARNING: trtime is '+str(state.aux_global['trtime'])+\
                ' but velocities reversed at time '+str(state.t)
    #Change to periodic BCs after initial pulse 
    if state.t>5*state.aux_global['tw1'] and solver.bc_lower[0]==0:
        solver.bc_lower[0]=2
        solver.bc_upper[0]=2


def zero_bc(grid,dim,t,qbc,mbc):
    """Set everything to zero"""
    print 'hello'
    if dim.nend==dim.n:
        qbc[:,-mbc:]=0.

def moving_wall_bc(grid,dim,t,qbc,mbc):
    """Initial pulse generated at left boundary by prescribed motion"""
    if dim.bc_lower==0:
        if dim.nstart==0:
           qbc[0,:mbc]=qbc[0,mbc] 
           t=state.t; t1=state.aux_global['t1']; tw1=state.aux_global['tw1']
           a1=state.aux_global['a1'];
           t0 = (t-t1)/tw1
           if abs(t0)<=1.: vwall = -a1*(1.+np.cos(t0*np.pi))
           else: vwall=0.
           for ibc in xrange(mbc-1):
               qbc[1,mbc-ibc-1] = 2*vwall*state.aux[1,ibc] - qbc[1,mbc+ibc]



def stegoton(use_petsc=0,kernel_language='Fortran',solver_type='classic',iplot=0,htmlplot=0,outdir='./_output'):
    """
    Stegoton problem.
    Nonlinear elasticity in periodic medium.
    See LeVeque & Yong (2003).

    $$\\epsilon_t - u_x = 0$$
    $$\\rho(x) u_t - \\sigma(\\epsilon,x)_x = 0$$
    """

    vary_Z=False

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()

    solver.kernel_language = kernel_language
    if kernel_language=='Python': solver.set_riemann_solver('nonlinear_elasticity')
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    #Use the same BCs for the aux array
    solver.aux_bc_lower = solver.bc_lower
    solver.aux_bc_upper = solver.bc_lower

    xlower=0.0; xupper=600.0
    cellsperlayer=6; mx=int(round(xupper-xlower))*cellsperlayer
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    grid = pyclaw.Grid(x)
    meqn = 2
    state = pyclaw.State(grid,meqn)

    #Set global parameters
    alpha = 0.5
    KA    = 1.0
    KB    = 4.0
    rhoA  = 1.0
    rhoB  = 4.0
    state.aux_global = {}
    state.aux_global['t1']    = 10.0
    state.aux_global['tw1']   = 10.0
    state.aux_global['a1']    = 0.0
    state.aux_global['alpha'] = alpha
    state.aux_global['KA'] = KA
    state.aux_global['KB'] = KB
    state.aux_global['rhoA'] = rhoA
    state.aux_global['rhoB'] = rhoB
    state.aux_global['trtime'] = 250.0
    state.aux_global['trdone'] = False

    #Initialize q and aux
    xc=grid.x.center
    state.aux=setaux(xc,rhoB,KB,rhoA,KA,alpha,solver.aux_bc_lower[0],xupper=xupper)
    qinit(state,ic=2,a2=1.0,xupper=xupper)

    tfinal=500.; nout = 10;

    solver.max_steps = 5000000
    solver.fwave = True 
    solver.start_step = b4step 
    solver.user_bc_lower=moving_wall_bc
    solver.user_bc_upper=zero_bc
    solver.mwaves=2

    if solver_type=='sharpclaw':
        solver.lim_type = 2
        solver.char_decomp=0

    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.outstyle = 1
    claw.nout = nout
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver


    if vary_Z==True:
        #Zvalues = [1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0]
        Zvalues = [3.5,4.0]
        #a2values= [0.9766161130681, 1.0888194560100042, 1.1601786315361329, 1.20973731651806, 1.2462158254919984]

        for ii,Z in enumerate(Zvalues):
            a2=1.0 #a2values[ii]
            KB = Z
            rhoB = Z
            state.aux_global['KB'] = KB
            state.aux_global['rhoB'] = rhoB
            state.aux_global['trdone'] = False
            state.aux=setaux(xc,rhoB,KB,rhoA,KA,alpha,bc_lower,xupper=xupper)
            grid.x.bc_lower=2
            grid.x.bc_upper=2
            state.t = 0.0
            qinit(state,ic=2,a2=a2)
            init_solution = Solution(state)
            claw.solution = init_solution
            claw.solution.t = 0.0

            claw.tfinal = tfinal
            claw.outdir = './_output_Z'+str(Z)+'_'+str(cellsperlayer)
            status = claw.run()

    else:
        # Solve
        status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(stegoton)
