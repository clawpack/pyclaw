#!/usr/bin/env python
# encoding: utf-8

from clawpack import pyclaw
import numpy as np
from clawpack.riemann import rp_advection
from clawpack.riemann import rp1_advection
import logging
from timeit import timeit

new_solvers = {}
try:
    import next.advection_module
    import next.classic1_module
    new_solvers['module'] = ('Fortran',
                             next.advection_module,
                             next.classic1_module)
except ImportError as err:
    print "Unable to import Module variant", err
    print "This module has not been implemented yet"
#    print "You may need to run python setup.py build_ext"

try:
    import next.advection_iso_c
    import next.classic1_iso_c
    new_solvers['iso_c'] = ('Fortran',
                            next.advection_iso_c,
                            next.classic1_iso_c)
except ImportError as err:
#    print "Unable to import ISO C variant", err
    print "This module has not been implemented yet"
    print "You may need to run python setup.py build_ext"

solvers = {
    'current_fortran' : ('Fortran', rp1_advection),
    'current_python'  : ('Python', rp_advection.rp_advection_1d)
}

solvers.update(new_solvers)

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
    riemann_compare.claw = claw


def verify(dx, q0, qfinal):
    if q0 is not None and qfinal is not None:
        test = dx*np.linalg.norm(qfinal-q0,1)
        return test
#        assert(abs(test-3.204e-04) < 1e-06, "Expected %e, got %e" % )

def compare(nx=1000):
    """
    Tests a variety of Riemann solver ideas on 1D advection
    """

    import riemann_compare

    solvers = riemann_compare.solvers

    times, tests = {}, {}

    for name, solver_tuple in solvers.iteritems():
        solver = pyclaw.ClawSolver1D()
        solver.kernel_language, solver.rp = solver_tuple[0:2]
        solver.num_waves = rp_advection.num_waves
        solver.bc_lower[0] = 2
        solver.bc_upper[0] = 2
        riemann_compare.solver = solver
        riemann_compare.state, riemann_compare.domain = init(nx)

        # benchmark
        t = timeit(stmt='riemann_compare.run()',
                   number=1, setup = 'import riemann_compare')
        times[name] = t

        # verify
        claw = riemann_compare.claw
        test = verify(claw.solution.domain.grid.delta[0],
                      claw.frames[0].state.get_q_global(),
                      claw.frames[claw.num_output_times].state.get_q_global())
        tests[name] = test

    return times, tests

if __name__=="__main__":
    import riemann_compare
    disable_loggers()

    nx = 500

    times, tests = compare(nx=nx)
    solvers = riemann_compare.solvers

    print "\nRiemann comparison on 1D advection with %d grid points" % nx

    print "\n=====TIME====="
    for name in solvers.keys():
        print "%-25s: %g" % (name, times[name])

    print "\n===ACCURACY==="
    for name in solvers.keys():
        print "%-25s: %g" % (name, tests[name])
