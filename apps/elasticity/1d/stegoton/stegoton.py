#!/usr/bin/env python
# encoding: utf-8
"""
stegoton.py

Stegoton problem.
Nonlinear elasticity in periodic medium.

$$\\epsilon_t - u_x = 0$$
$$\\rho(x) u_t - \\sigma(\\epsilon,x)_x = 0$$
"""


import os
import numpy as np

from petclaw.grid import PCDimension as Dimension
from petclaw.grid import PCGrid as Grid
from pyclaw.solution import Solution
from pyclaw.controller import Controller
from petclaw.evolve.petclaw import PetClawSolver1D
from petsc4py import PETSc

def qinit(grid):

    # Initial Data parameters
    a1 = grid.aux_global['a1']

    x =grid.x.center
    q=np.zeros([len(x),grid.meqn])
    
    # Gaussian
    mbc=grid.mbc
    sigma = a1*np.exp(-((x-75.)/10.)**2.)
    q[:,0] = np.log(sigma+1.)/grid.aux[mbc:-mbc,1]
    q[:,0] = 0.

    grid.q=q


def setaux(x,rhoB=4,KB=4,rhoA=1,KA=1,alpha=0.5,xlower=0.,xupper=150.):
    aux=np.empty([len(x),2])

    xfrac = x-np.floor(x)

    #Density:
    aux[:,0] = rhoA*(xfrac<alpha)+rhoB*(xfrac>=alpha)
    #Bulk modulus:
    aux[:,1] = KA  *(xfrac<alpha)+KB  *(xfrac>=alpha)
    for i,xi in enumerate(x):
        if xi<xlower:
            aux[i,0]=rhoA
            aux[i,1]=KA
        if xi>xupper:
            aux[i,0]=rhoB
            aux[i,1]=KB
    return aux

    
def b4step(solver,solutions):
    #Reverse velocity at trtime
    #Note that trtime should be an output point
    grid = solutions['n'].grids[0]
    if grid.t>=grid.aux_global['trtime']-1.e-10 and not grid.aux_global['trdone']:
        print 'Time reversing'
        grid.q[:,1]=-grid.q[:,1]
        solutions['n'].grid.q=grid.q
        grid.aux_global['trdone']=True
        if grid.t>grid.aux_global['trtime']:
            print 'WARNING: trtime is '+str(grid.aux_global['trtime'])+\
                ' but velocities reversed at time '+str(grid.t)
    #Change to periodic BCs after initial pulse 
    if grid.t>5*grid.aux_global['tw1'] and grid.x.mthbc_lower==0:
        grid.x.mthbc_lower=2
        grid.x.mthbc_upper=2
        #Also make aux arrays periodic
        xghost=grid.x.centerghost
        for i,x in enumerate(xghost):
            if x<grid.x.lower:
                grid.aux[i,0]=grid.aux_global['rhoB']
                grid.aux[i,1]=grid.aux_global['KB']
            if x>grid.x.upper:
                grid.aux[i,0]=grid.aux_global['rhoA']
                grid.aux[i,1]=grid.aux_global['KA']


    
class PPCGrid(Grid):
    def init_q_petsc_structures(self):
        r"""
        Override usual behavior in order to have BCs that start as non-periodic
        and then switch to periodic.
        """
        periodic_type = PETSc.DA.PeriodicType.X

        self.q_da = PETSc.DA().create(dim=self.ndim,
                                    dof=self.meqn,
                                    sizes=self.n, 
                                    periodic_type = periodic_type,
                                    stencil_width=self.mbc,
                                    comm=PETSc.COMM_WORLD)
        self.gqVec = self.q_da.createGlobalVector()
        self.lqVec = self.q_da.createLocalVector()

        #Now set up the local indices:
        ranges = self.q_da.getRanges()
        for i,range in enumerate(ranges):
            self.dimensions[i].nstart=range[0]
            self.dimensions[i].nend  =range[1]
            
 
def zero_bc(grid,dim,qbc):
    """Set everything to zero"""
    if dim.mthbc_upper==0:
        if dim.centerghost[-1]>150.:
           qbc[-grid.mbc:,:]=0.

def moving_wall_bc(grid,dim,qbc):
    """Initial pulse generated at left boundary by prescribed motion"""
    if dim.mthbc_lower==0:
        if dim.centerghost[0]<0:
           qbc[:grid.mbc,0]=qbc[grid.mbc,0] 
           t=grid.t; t1=grid.aux_global['t1']; tw1=grid.aux_global['tw1']
           a1=grid.aux_global['a1']; mbc=grid.mbc
           t0 = (t-t1)/tw1
           if abs(t0)<=1.: vwall = -a1*(1.+np.cos(t0*np.pi))
           else: vwall=0.
           for ibc in xrange(mbc-1):
               qbc[mbc-ibc-1,1] = 2*vwall*grid.aux[ibc,1] - qbc[mbc+ibc,1]


if __name__ == "__main__":
    import time
    start=time.time()
    # Initialize grids and solutions
    xlower=0.0; xupper=150.0
    cellsperlayer=24; mx=150*cellsperlayer
    x = Dimension('x',xlower,xupper,mx,mthbc_lower=0,mthbc_upper=0,mbc=2)
    grid = PPCGrid(x)
    grid.meqn = 2
    grid.t = 0.0

    #Set global parameters
    alpha = 0.5
    KA    = 1.0
    KB    = 4.0
    rhoA  = 1.0
    rhoB  = 4.0
    grid.aux_global = {}
    grid.aux_global['t1']    = 10.0
    grid.aux_global['tw1']   = 10.0
    grid.aux_global['a1']    = 0.1
    grid.aux_global['alpha'] = alpha
    grid.aux_global['KA'] = KA
    grid.aux_global['KB'] = KB
    grid.aux_global['rhoA'] = rhoA
    grid.aux_global['rhoB'] = rhoB
    grid.aux_global['trtime'] = 250.0
    grid.aux_global['trdone'] = False

    # Initilize petsc Structures
    grid.init_q_petsc_structures()
        
    #Initialize q and aux
    xghost=grid.x.centerghost
    grid.aux=setaux(xghost,rhoB,KB,rhoA,KA,alpha)
    qinit(grid)
    init_solution = Solution(grid)

    Kmax=max(grid.aux_global['KA'],grid.aux_global['KB'])
    emax=np.max(grid.q[:,0])
    smax=np.sqrt(Kmax*np.exp(Kmax*emax)) #This isn't quite right

    # Solver setup
    solver = PetClawSolver1D(kernelsType = 'F')

    tfinal=500.; nout = 10; tout=tfinal/nout
    dt_rough = 1.45*grid.x.d/smax
    nsteps = np.ceil(tout/dt_rough)
    solver.dt = tout/nsteps

    solver.max_steps = 5000
    solver.set_riemann_solver('nel')
    solver.order = 2
    solver.mthlim = [4,4]
    solver.dt_variable = False
    solver.fwave = True 
    solver.start_step = b4step 
    solver.user_bc_lower=moving_wall_bc
    solver.user_bc_upper=zero_bc

    use_controller = True

    if(use_controller):

    # Controller instantiation
        claw = Controller()
        claw.outdir = './_output'
        claw.keep_copy = False
        claw.nout = nout
        claw.outstyle = 1
        claw.output_format = 'petsc'
        claw.tfinal = tfinal
        claw.solutions['n'] = init_solution
        claw.solver = solver

        # Solve
        status = claw.run()
        end=time.time()
        print 'job took '+str(end-start)+' seconds'


    else:
        sol = {"n":init_solution}
        
        solver.evolve_to_time(sol,.4)
        sol = sol["n"]
