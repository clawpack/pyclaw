#!/usr/bin/env python
# encoding: utf-8

from clawpack import pyclaw
import numpy as np
from clawpack.riemann import rp_advection
from clawpack.riemann import rp1_advection
import logging
from timeit import timeit

def disable_loggers():
    root_logger = logging.getLogger()
    root_logger.disabled = True

def init(nx):
    x = pyclaw.Dimension('x',0.0,1.0,nx)
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)
    state.problem_data['u']=1.

    grid = state.grid
    xc=grid.x.centers
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    return state, domain

def run():
    """
    Given riemann_compare module variables state, domain, solver, runs an
    instance of pyclaw.Controller until tfinal=1
    """
    import riemann_compare

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(riemann_compare.state,
                                    riemann_compare.domain)
    claw.solver = riemann_compare.solver

    claw.output_format = None

    claw.tfinal =1.0
    status = claw.run()
    return claw

def compare(nx=1000):
    """
    Tests a variety of Riemann solver ideas on 1D advection
    """

    import riemann_compare

    solvers = {
        'current_python'  : ('Python', rp_advection.rp_advection_1d),
        'current_fortran' : ('Fortran', rp1_advection),
    }
    times = {}

    for name, solver_tuple in solvers.iteritems():
        solver = pyclaw.ClawSolver1D()
        solver.kernel_language, solver.rp = solver_tuple
        solver.num_waves = rp_advection.num_waves
        solver.bc_lower[0] = 2
        solver.bc_upper[0] = 2
        riemann_compare.solver = solver
        riemann_compare.state, riemann_compare.domain = init(nx)

        t = timeit(stmt='riemann_compare.run()',
                   number=1, setup = 'import riemann_compare')
        times[name] = t
    return solvers, times

if __name__=="__main__":
    disable_loggers()
    solvers, times = compare()

    for name, solver_tuple in solvers.iteritems():
        print "%-25s: %g" % (name, times[name])
