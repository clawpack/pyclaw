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

    grid.q=q


def setaux(grid):
    maux = 2
    xghost=grid.x.centerghost
    grid.aux=np.empty([len(xghost),maux])

    alpha=grid.aux_global['alpha']
    KA = grid.aux_global['KA']
    KB = grid.aux_global['KB']
    rhoA = grid.aux_global['rhoA']
    rhoB = grid.aux_global['rhoB']

    xfrac = xghost-np.floor(xghost)

    #Density:
    grid.aux[:,0] = rhoA*(xfrac<alpha)+rhoB*(xfrac>=alpha)
    #Bulk modulus:
    grid.aux[:,1] = KA  *(xfrac<alpha)+KB  *(xfrac>=alpha)
    for i,x in enumerate(xghost):
        if x<grid.x.lower:
            grid.aux[i,0]=rhoA
            grid.aux[i,1]=KA
        if x>grid.x.upper:
            grid.aux[i,0]=rhoB
            grid.aux[i,1]=KB
    
    
def b4step(solver,solutions):
    #Reverse velocity at trtime
    #Note that trtime should be an output point
    grid = solutions['n'].grids[0]
    if grid.t>=grid.aux_global['trtime'] and not grid.aux_global['trdone']:
        grid.q[:,1]=-grid.q[:,1]    
        if grid.t>grid.aux_global['trtime']:
            print 'WARNING: trtime is '+str(grid.aux_global['trtime'])+\
                ' but velocities reversed at time '+str(grid.t)
    #Change to periodic BCs after initial pulse 
    if grid.t>5*grid.aux_global['tw1'] and grid.x.mthbc_lower==0:
        grid.x.mthbc_lower=2
        grid.x.mthbc_upper=2
    

# Initialize grids and solutions
xlower=0.0; xupper=150.0
cellsperlayer=24; mx=150*cellsperlayer
x = Dimension('x',xlower,xupper,mx,mthbc_lower=1,mthbc_upper=1,mbc=2)
grid = Grid(x)
grid.meqn = 2
grid.t = 0.0

#Set global parameters
grid.aux_global = {}
grid.aux_global['t1']    = 10.0
grid.aux_global['tw1']   = 10.0
grid.aux_global['a1']    = 0.2
grid.aux_global['alpha'] = 0.5
grid.aux_global['KA']    = 1.0
grid.aux_global['KB']    = 4.0
grid.aux_global['rhoA']  = 1.0
grid.aux_global['rhoB']  = 4.0
grid.aux_global['trtime'] = 600.0
grid.aux_global['trdone'] = False

# Initilize petsc Structures
grid.init_q_petsc_structures()
    
#Initialize q and aux
setaux(grid)
qinit(grid)
init_solution = Solution(grid)

Kmax=max(grid.aux_global['KA'],grid.aux_global['KB'])
emax=np.max(grid.q[:,0])
smax=Kmax*np.exp(Kmax*emax)

# Solver setup
solver = PetClawSolver1D(kernelsType = 'P')
solver.dt = 0.2*grid.x.d/smax
print solver.dt
solver.max_steps = 5000
solver.set_riemann_solver('nel')
solver.order = 2
solver.mthlim = [4,4]
solver.dt_variable = False #Amal: need to handle the case dt_variable.
solver.fwave = True 
solver.start_step = b4step 

use_controller = True

if(use_controller):

# Controller instantiation
    claw = Controller()
    claw.outdir = './_output'
    claw.keep_copy = False
    claw.nout = 100
    claw.outstyle = 1
    claw.output_format = 'petsc'
    claw.tfinal = 250.0
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()


else:
    sol = {"n":init_solution}
    
    solver.evolve_to_time(sol,.4)
    sol = sol["n"]
