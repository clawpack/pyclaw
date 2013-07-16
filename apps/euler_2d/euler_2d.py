#!/usr/bin/env python
# encoding: utf-8

"""
Solve the Euler equations of compressible fluid dynamics.
"""
from clawpack import pyclaw
from clawpack import riemann

solver = pyclaw.ClawSolver2D(riemann.rp2_euler_4wave)
solver.all_bcs = pyclaw.BC.extrap

domain = pyclaw.Domain([0.,0.],[1.,1.],[100,100])
solution = pyclaw.Solution(solver.num_eqn,domain)
gamma = 1.4
solution.problem_data['gamma']  = gamma

# Set initial data
xx,yy = domain.grid.p_centers
l = xx<0.5; r = xx>=0.5; b = yy<0.5; t = yy>=0.5
solution.q[0,...] = 2.*l*t + 1.*l*b + 1.*r*t + 3.*r*b
solution.q[1,...] = 0.75*t - 0.75*b
solution.q[2,...] = 0.5*l  - 0.5*r
solution.q[3,...] = 0.5*solution.q[0,...]*(solution.q[1,...]**2+solution.q[2,...]**2) + 1./(gamma-1.)

#solver.evolve_to_time(solution,tend=0.3)
claw = pyclaw.Controller()
claw.tfinal = 0.3
claw.solution = solution
claw.solver = solver

status = claw.run()

#pyclaw.plot.interactive_plot()
