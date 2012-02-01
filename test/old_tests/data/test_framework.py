#!/usr/bin/env python

import pyclaw
import numpy as np
        
def true_solution(t):
    # Setup problem parameters
    beta = 100.0
    gamma = 0.0
    x0 = 0.3
    x1 = 0.6
    x2 = 0.8
    efix = False

    # Create empty patch
    x = pyclaw.solution.Dimension('x',0.0,1.0,100)
    patch = pyclaw.solution.Patch(x)
    patch.t = t
    patch.problem_data['efix'] = efix
    patch.empty_q()

    # Gaussian
    qg = np.exp(-beta * (x.centers-x0)**2) * np.cos(gamma * (x.centers - x0))

    # Step Function
    qs = (x.centers > x1) * 1.0 - (x.centers > x2) * 1.0


    patch.q[0,:] = qg + qs
    
    return pyclaw.solution.Solution(patch)
        
# Default settings
controller = pyclaw.controller.Controller()
controller.verbosity = 0
controller.output_format = None
controller.keep_copy = True
controller.output_style = 2
controller.out_times = [0.0,0.5,1.0]
        
# Setup solver, using mainly defaults so should be setup later
controller.solver = pyclaw.evolve.clawpack.ClawSolver1D()
controller.solver.set_riemann_solver('burgers')
controller.solver.mthlim = 3
controller.solver.order = 2
controller.solver.dt = 0.0001
controller.solver.max_steps = 5000
controller.solver.num_waves = 1

controller.solution = true_solution(0.0)
controller.run()

controller.solution.write(0,path='./burgers_test/')

import matplotlib.pyplot as plt
for (i,sol) in enumerate(controller.frames):
    plt.figure(1)
    plt.subplot(3,2,i+1)
    plt.plot(sol.patchs[0].dimensions[0].centers,sol.patchs[0].q[0,:],'k')
    if not(i == 0):
        sol.write(i,path="./burgers_test")
plt.show()

# for (i,solution) in enumerate(sol):
#     solution.write(i,path='./limiter_test')
