#!/usr/bin/env python
# encoding: utf-8

from clawpack import pyclaw
import numpy as np
import logging

def debug_loggers():
    """
    Turn on maximum debugging from all loggers.
    """

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_logger.propagate = True


def disable_loggers():
    """
    Disable all loggers (quiet runs)
    """

    root_logger = logging.getLogger()
    root_logger.disabled = True

from clawpack.riemann import advection_1D
fsolver_1D = pyclaw.ClawSolver1D(advection_1D)
fsolver_1D.kernel_language = 'Fortran'

from clawpack.riemann import advection_1D_py
pysolver_1D = pyclaw.ClawSolver1D(advection_1D_py.advection_1D)
pysolver_1D.kernel_language = 'Python'

from clawpack.riemann import shallow_roe_with_efix_2D
fsolver_2D = pyclaw.ClawSolver2D(shallow_roe_with_efix_2D)
fsolver_2D.kernel_language = 'Fortran'

solvers_1D = {
    'current_fortran' : fsolver_1D,
    'current_python'  : pysolver_1D
}

solvers_2D = {
    'current_fortran' : fsolver_2D
}

new_solvers_1D = {}

# Here we try and bring in "experimental" solvers for comparison
# If we can't bring in the solver, complain and move on...

try:
    from clawpack.pyclaw.examples import iso_c_advection

    iso_c_solver = iso_c_advection.ISO_C_ClawSolver1D(None, 'clawpack.pyclaw')

    iso_c_solver.kernel_language = 'Fortran'
    iso_c_solver.rp = iso_c_advection.iso_c_rp1_advection(1.0)
    iso_c_solver.num_waves = 1

    new_solvers_1D['iso_c'] = iso_c_solver

except (ImportError, OSError) as err:
    print("Unable to import ISO C variant", err)

solvers_1D.update(new_solvers_1D)


def verify_1D(dx, q0, qfinal):
    if q0 is not None and qfinal is not None:
        test = dx*np.linalg.norm(qfinal-q0,1)
        return test

def verify_2D(dx, qfinal):
    if qfinal is not None:
        test = dx*np.linalg.norm(qfinal)
        return test

def compare_1D(nx=1000):
    """
    Tests a variety of Riemann solver ideas on 1D advection
    """

    from . import compare_solvers
    import time

    solvers = compare_solvers.solvers_1D

    times, tests = {}, {}

    for name, solver in solvers.items():
        solver.bc_lower[0] = pyclaw.BC.periodic
        solver.bc_upper[0] = pyclaw.BC.periodic

        from clawpack.pyclaw.examples.advection_1d import advection
        claw = advection.advection(nx=nx,outdir=compare_solvers.outdir)
        claw.solver = solver
        claw.keep_copy = True

        # benchmark
        t0 = time.clock()
        claw.run()
        t1 = time.clock()
        t = t1-t0
        times[name] = t

        # verify
        test = verify_1D(claw.solution.domain.grid.delta[0],
                         claw.frames[0].state.get_q_global(),
                         claw.frames[claw.num_output_times].state.get_q_global())
        tests[name] = test

    return times, tests


def compare_2D(nx=(250,250)):
    """
    Tests a variety of Riemann solver ideas on 2D shallow water equation
    """
    from . import compare_solvers
    import time

    solvers = compare_solvers.solvers_2D

    times, tests = {}, {}

    for name, solver in solvers.items():
        solver.num_waves = 3
        solver.bc_lower[0] = pyclaw.BC.extrap
        solver.bc_upper[0] = pyclaw.BC.wall
        solver.bc_lower[1] = pyclaw.BC.extrap
        solver.bc_upper[1] = pyclaw.BC.wall

        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=1

        from clawpack.pyclaw.examples.shallow_2d import shallow2D
        claw = shallow2D.shallow2D(outdir=compare_solvers.outdir)
        claw.solver = solver
        claw.keep_copy = True

        t0 = time.clock()
        claw.run()
        t1 = time.clock()
        t = t1-t0
        times[name] = t

        # verify
        test = verify_2D(claw.solution.domain.grid.delta[0],
                         claw.frames[claw.num_output_times].state.get_q_global())
        tests[name] = test

    return times, tests

if __name__=="__main__":
    from . import compare_solvers
    disable_loggers()
#    debug_loggers()

    import sys

    vis = True
    
    if vis:
        compare_solvers.outdir='./_output'
    else:
        compare_solvers.outdir = None
        
    if len(sys.argv) > 1:
        nx_1D = int(sys.argv[1])
    else:
        nx_1D = 500

    if len(sys.argv) > 2:
        # '(2,2)' -> (2,2)
        nx_2D = tuple(int(i) for i in argv.split(','))
    else:
        nx_2D = 100,100

    def print_time_accuracy(times, tests, solvers):
        print("\n=====TIME=====")
        for name in solvers.keys():
            print("%-25s: %g" % (name, times[name]))

        print("\n===ACCURACY===")
        for name in solvers.keys():
            print("%-25s: %g" % (name, tests[name]))

    compare_solvers.tfinal = 1.0
    print("\nRiemann comparison on 1D advection to t=%g with %d grid points" % \
        (compare_solvers.tfinal, nx_1D))

    times, tests = compare_1D(nx=nx_1D)

    print_time_accuracy(times, tests, compare_solvers.solvers_1D)

    compare_solvers.tfinal = 2.5
    print(("\nRiemann comparison on 2D shallow water equation to t=%g" + \
               " with %dx%d grid points") % ((compare_solvers.tfinal,) + nx_2D))

    times, tests = compare_2D(nx=nx_2D)

    print_time_accuracy(times, tests, compare_solvers.solvers_2D)
