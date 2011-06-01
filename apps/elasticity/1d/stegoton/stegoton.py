#!/usr/bin/env python
# encoding: utf-8
import numpy as np

def qinit(grid,ic=2,a2=1.0,xupper=600.):
    x =grid.x.center
    grid.zeros_q()
    
    if ic==1: #Zero ic
        pass
    elif ic==2:
        # Gaussian
        mbc=grid.mbc
        sigma = a2*np.exp(-((x-xupper/2.)/10.)**2.)
        grid.q[0,:] = np.log(sigma+1.)/grid.aux[1,mbc:-mbc]


def setaux(x,rhoB=4,KB=4,rhoA=1,KA=1,alpha=0.5,xlower=0.,xupper=600.,mthbc=2):
    aux=np.zeros([3,len(x)],order='F')

    xfrac = x-np.floor(x)

    #Density:
    aux[0,:] = rhoA*(xfrac<alpha)+rhoB*(xfrac>=alpha)
    #Bulk modulus:
    aux[1,:] = KA  *(xfrac<alpha)+KB  *(xfrac>=alpha)
    for i,xi in enumerate(x):
        if xi<xlower:
            if mthbc==2:
                aux[0,i]=rhoB
                aux[1,i]=KB
            else:
                aux[0,i]=rhoA
                aux[1,i]=KA
        if xi>xupper:
            if mthbc==2:
                aux[0,i]=rhoA
                aux[1,i]=KA
            else:
                aux[0,i]=rhoB
                aux[1,i]=KB
    return aux

    
def b4step(solver,solutions):
    #Reverse velocity at trtime
    #Note that trtime should be an output point
    grid = solutions['n'].grids[0]
    if grid.t>=grid.aux_global['trtime']-1.e-10 and not grid.aux_global['trdone']:
        #print 'Time reversing'
        grid.q[1,:]=-grid.q[1,:]
        solutions['n'].grid.q=grid.q
        grid.aux_global['trdone']=True
        if grid.t>grid.aux_global['trtime']:
            print 'WARNING: trtime is '+str(grid.aux_global['trtime'])+\
                ' but velocities reversed at time '+str(grid.t)
    #Change to periodic BCs after initial pulse 
    if grid.t>5*grid.aux_global['tw1'] and solver.mthbc_lower[0]==0:
        solver.mthbc_lower[0]=2
        solver.mthbc_upper[0]=2
        #Also make aux arrays periodic
        xghost=grid.x.centerghost
        for i,x in enumerate(xghost):
            if x<grid.x.lower:
                grid.aux[0,i]=grid.aux_global['rhoB']
                grid.aux[1,i]=grid.aux_global['KB']
            if x>grid.x.upper:
                grid.aux[0,i]=grid.aux_global['rhoA']
                grid.aux[1,i]=grid.aux_global['KA']


def zero_bc(grid,dim,t,qbc):
    """Set everything to zero"""
    if dim.centerghost[-1]>150.:
        qbc[:,-grid.mbc:]=0.

def moving_wall_bc(grid,dim,t,qbc):
    """Initial pulse generated at left boundary by prescribed motion"""
    if dim.mthbc_lower==0:
        if dim.centerghost[0]<0:
           qbc[0,:grid.mbc]=qbc[0,grid.mbc] 
           t=grid.t; t1=grid.aux_global['t1']; tw1=grid.aux_global['tw1']
           a1=grid.aux_global['a1']; mbc=grid.mbc
           t0 = (t-t1)/tw1
           if abs(t0)<=1.: vwall = -a1*(1.+np.cos(t0*np.pi))
           else: vwall=0.
           for ibc in xrange(mbc-1):
               qbc[1,mbc-ibc-1] = 2*vwall*grid.aux[1,ibc] - qbc[1,mbc+ibc]



def stegoton(use_petsc=True,kernel_language='Fortran',solver_type='classic',iplot=False,htmlplot=False,outdir='./_output'):
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
    if kernel_language=='Python': solver.set_riemann_solver('nel')
    solver.mthbc_lower[0] = 2
    solver.mthbc_upper[0] = 2

    xlower=0.0; xupper=600.0
    cellsperlayer=6; mx=int(round(xupper-xlower))*cellsperlayer
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    grid = pyclaw.Grid(x)
    grid.meqn = 2
    grid.mbc=solver.mbc

    #Set global parameters
    alpha = 0.5
    KA    = 1.0
    KB    = 4.0
    rhoA  = 1.0
    rhoB  = 4.0
    grid.aux_global = {}
    grid.aux_global['t1']    = 10.0
    grid.aux_global['tw1']   = 10.0
    grid.aux_global['a1']    = 0.0
    grid.aux_global['alpha'] = alpha
    grid.aux_global['KA'] = KA
    grid.aux_global['KB'] = KB
    grid.aux_global['rhoA'] = rhoA
    grid.aux_global['rhoB'] = rhoB
    grid.aux_global['trtime'] = 250.0
    grid.aux_global['trdone'] = False

    #Initialize q and aux
    xghost=grid.x.centerghost
    grid.aux=setaux(xghost,rhoB,KB,rhoA,KA,alpha,solver.mthbc_lower[0],xupper=xupper)
    qinit(grid,ic=2,a2=1.0,xupper=xupper)

    Kmax=max(grid.aux_global['KA'],grid.aux_global['KB'])
    emax=2*grid.aux_global['a1']*np.sqrt(Kmax)
    smax=np.sqrt(np.exp(Kmax*emax)) #Works only for K=rho

    tfinal=500.; nout = 10; tout=tfinal/nout
    dt_rough = 0.5*grid.x.d/smax 
    nsteps = np.ceil(tout/dt_rough)
    solver.dt = tout/nsteps/2.

    solver.max_steps = 5000000
    solver.fwave = True 
    solver.start_step = b4step 
    solver.user_bc_lower=moving_wall_bc
    solver.user_bc_upper=zero_bc
    solver.mwaves=2

    if solver_type=='sharpclaw':
        solver.lim_type = 2
        solver.time_integrator='SSP33'
        solver.char_decomp=0

    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.outstyle = 1
    claw.nout = nout
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(grid)
    claw.solver = solver


    if vary_Z==True:
        #Zvalues = [1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0]
        Zvalues = [3.5,4.0]
        #a2values= [0.9766161130681, 1.0888194560100042, 1.1601786315361329, 1.20973731651806, 1.2462158254919984]

        for ii,Z in enumerate(Zvalues):
            a2=1.0 #a2values[ii]
            KB = Z
            rhoB = Z
            grid.aux_global['KB'] = KB
            grid.aux_global['rhoB'] = rhoB
            grid.aux_global['trdone'] = False
            grid.aux=setaux(xghost,rhoB,KB,rhoA,KA,alpha,mthbc_lower,xupper=xupper)
            grid.x.mthbc_lower=2
            grid.x.mthbc_upper=2
            grid.t = 0.0
            qinit(grid,ic=2,a2=a2)
            init_solution = Solution(grid)
            claw.solutions['n'] = init_solution
            claw.solutions['n'].t = 0.0

            claw.tfinal = tfinal
            claw.outdir = './_output_Z'+str(Z)+'_'+str(cellsperlayer)
            status = claw.run()

    else:
        # Solve
        status = claw.run()

    if htmlplot:  pyclaw.plot.plotHTML(outdir=outdir)
    if iplot:     pyclaw.plot.plotInteractive(outdir=outdir)


if __name__=="__main__":
    import sys
    from pyclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    stegoton(*args,**kwargs)
