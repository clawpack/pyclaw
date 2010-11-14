#!/usr/bin/python
#!/usr/bin/env python
# encoding: utf-8
"""
advection.py

Example python script for solving the 1d advection equation.
"""

import os, sys



try:
    import numpy as np
    from petsc4py import PETSc
    
except:
    sys.path.append("/scratch/ketch/ksl/petsc4py/dev-aug29/ppc450d/lib/python/")
    sys.path.append("/scratch/ketch/ksl/numpy/dev-aug29/ppc450d/lib/python/")
    
    import numpy as np
    from petsc4py import PETSc

try:
    from petclaw.grid import PCDimension as Dimension
    from petclaw.grid import PCGrid as Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.petclaw import PetClawSolver1D
    from pyclaw.controller import Controller
except:
    #PATH = os.environ['PATH']
    #print "Path::",PATH

    #CLAW = os.environ['CLAW']
    #PETCLAW = os.environ['PETCLAW']
    #PYCLAW = CLAW+"/python"
    #PETCLAW_PKG = PETCLAW + "/src"
    #sys.path.append(PETCLAW_PKG)
    #sys.path.append(PYCLAW)
    sys.path.append("/scratch/ketch/petclaw/src/")
    sys.path.append("/scratch/ketch/clawpack4petclaw/python/")
    from petclaw.grid import PCDimension as Dimension
    from petclaw.grid import PCGrid as Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.petclaw import PetClawSolver1D
    from pyclaw.controller import Controller



def qinit(grid):

    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Initial Data parameters
    ic = 3
    beta = 100.
    gamma = 0.
    x0 = 0.3
    x1 = 0.7
    x2 = 0.9
    
    # Create an array with fortran native ordering
    x =grid.x.center
   
    q=np.zeros([len(x),grid.meqn])
    
    # Gaussian
    qg = np.exp(-beta * (x-x0)**2) * np.cos(gamma * (x - x0))

    # Step Function
    qs = (x > x1) * 1.0 - (x > x2) * 1.0
    
    if ic == 1: q[:,0] = qg
    elif ic == 2: q[:,0] = qs
    elif ic == 3: q[:,0] = qg + qs
    grid.q=q


# Initialize grids and solutions
numprocs=1
mx=4096*numprocs
x = Dimension('x',0.0,1.0,mx,mthbc_lower=2,mthbc_upper=2)
grid = Grid(x)
grid.aux_global['u']=-1.
grid.meqn = 1
grid.t = 0.0
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = PetClawSolver1D(kernelsType = 'F')

solver.dt = 0.8*grid.x.d
solver.max_steps = 2000

solver.set_riemann_solver('advection')
solver.order = 2
solver.mthlim = 4
solver.dt_variable = False

tfinal = 500*solver.dt

useController = False


if useController:

    # Controller instantiation
    claw = Controller()
    claw.outdir = './_output/'
    claw.keep_copy = False
    claw.nout = 1
    claw.outstyle = 1
    claw.output_format = 'petsc'
    claw.tfinal = tfinal
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()

else:
    import time
    sol = {"n":init_solution}
    start=time.time()
    solver.evolve_to_time(sol,tfinal)
    end=time.time()
    duration1 = end-start
    print 'First part of job took '+str(duration1)+' seconds'

    start=time.time()
    solver.evolve_to_time(sol,tfinal*2)
    end=time.time()
    duration2 = end-start
    print 'Second part of job took '+str(duration2)+' seconds'
