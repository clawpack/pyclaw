#!/usr/bin/env python
# encoding: utf-8

from clawpack import pyclaw
import numpy as np
from clawpack.riemann import rp_advection
from clawpack.riemann import rp1_advection
from clawpack.riemann import rp2_shallow_roe_with_efix
import logging
from timeit import timeit

def debug_loggers():
    """
    Turn on maximimum debugging from all loggers.
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

fsolver_1D = pyclaw.ClawSolver1D()
fsolver_1D.kernel_language = 'Fortran'
fsolver_1D.rp = rp1_advection

fsolver_2D = pyclaw.ClawSolver2D()
fsolver_2D.kernel_language = 'Fortran'
fsolver_2D.rp = rp2_shallow_roe_with_efix

pysolver_1D = pyclaw.ClawSolver1D()
pysolver_1D.kernel_language = 'Python'
pysolver_1D.rp = rp_advection.rp_advection_1d

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
    import next

    iso_c_solver = next.ISO_C_ClawSolver1D(None, 'clawpack.pyclaw')

    iso_c_solver.kernel_language = 'Fortran'
    iso_c_solver.rp = next.iso_c_rp1_advection(1.0)

    new_solvers_1D['iso_c'] = iso_c_solver

except (ImportError, OSError) as err:
    print "Unable to import ISO C variant", err

solvers_1D.update(new_solvers_1D)

def init_1D(nx):
    """ Initialize the grid for the 1D advection problem"""
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

def init_2D(nx):
    """ Initialize the grid for the 2D shallow water problem"""

    # Domain:
    xlower = -2.5
    xupper = 2.5
    mx = nx[0]
    ylower = -2.5
    yupper = 2.5
    my = nx[1]
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    domain = pyclaw.Domain([x,y])

    num_eqn = 3  # Number of equations
    state = pyclaw.State(domain,num_eqn)

    grav = 1.0 # Parameter (global auxiliary variable)
    state.problem_data['grav'] = grav

    # Initial solution
    # ================
    # Riemann states of the dam break problem
    radDam = 0.5
    hl = 2.
    ul = 0.
    vl = 0.
    hr = 1.
    ur = 0.
    vr = 0.
    x0=0.
    y0=0.
    xCenter = state.grid.x.centers
    yCenter = state.grid.y.centers
    Y,X = np.meshgrid(yCenter,xCenter)
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    state.q[0,:,:] = hl*(r<=radDam) + hr*(r>radDam)
    state.q[1,:,:] = hl*ul*(r<=radDam) + hr*ur*(r>radDam)
    state.q[2,:,:] = hl*vl*(r<=radDam) + hr*vr*(r>radDam)

    return state, domain

def run():
    """
    Given riemann_compare module variables state, domain, solver, runs an
    instance of pyclaw.Controller until riemann_compare.tfinal
    """

    import riemann_compare

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(riemann_compare.state,
                                    riemann_compare.domain)
    claw.solver = riemann_compare.solver

    claw.output_format = None

    claw.tfinal = riemann_compare.tfinal
    status = claw.run()
    riemann_compare.claw = claw


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

    import riemann_compare

    solvers = riemann_compare.solvers_1D

    times, tests = {}, {}

    for name, solver in solvers.iteritems():
        solver.num_waves = rp_advection.num_waves
        solver.bc_lower[0] = 2
        solver.bc_upper[0] = 2
        riemann_compare.solver = solver
        riemann_compare.state, riemann_compare.domain = init_1D(nx)

        # benchmark
        t = timeit(stmt='riemann_compare.run()',
                   number=1, setup = 'import riemann_compare')
        times[name] = t

        # verify
        claw = riemann_compare.claw
        test = verify_1D(claw.solution.domain.grid.delta[0],
                         claw.frames[0].state.get_q_global(),
                         claw.frames[claw.num_output_times].state.get_q_global())
        tests[name] = test

    return times, tests


def compare_2D(nx=(250,250)):
    """
    Tests a variety of Riemann solver ideas on 2D shallow water equation
    """

    import riemann_compare

    solvers = riemann_compare.solvers_2D

    times, tests = {}, {}

    for name, solver in solvers.iteritems():

        solver.num_waves = 3
        solver.bc_lower[0] = pyclaw.BC.extrap
        solver.bc_upper[0] = pyclaw.BC.wall
        solver.bc_lower[1] = pyclaw.BC.extrap
        solver.bc_upper[1] = pyclaw.BC.wall

        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=1

#        solver.dt_initial = 0.001
#        solver.dt_variable = False

        riemann_compare.solver = solver
        riemann_compare.state, riemann_compare.domain = init_2D(nx)

        # benchmark
        t = timeit(stmt='riemann_compare.run()',
                   number=1, setup = 'import riemann_compare')
        times[name] = t

        # verify
        claw = riemann_compare.claw
        test = verify_2D(claw.solution.domain.grid.delta[0],
                         claw.frames[claw.num_output_times].state.get_q_global())
        tests[name] = test

    return times, tests

if __name__=="__main__":
    import riemann_compare
    disable_loggers()
#    debug_loggers()

    import sys

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
        print "\n=====TIME====="
        for name in solvers.keys():
            print "%-25s: %g" % (name, times[name])

        print "\n===ACCURACY==="
        for name in solvers.keys():
            print "%-25s: %g" % (name, tests[name])

    riemann_compare.tfinal = 1.0
    print "\nRiemann comparison on 1D advection to t=%g with %d grid points" % \
        (riemann_compare.tfinal, nx_1D)

    times, tests = compare_1D(nx=nx_1D)

    print_time_accuracy(times, tests, riemann_compare.solvers_1D)

    riemann_compare.tfinal = 2.5
    print ("\nRiemann comparison on 2D shallow water equation to t=%g" + \
               " with %dx%d grid points") % ((riemann_compare.tfinal,) + nx_2D)

    times, tests = compare_2D(nx=nx_2D)

    print_time_accuracy(times, tests, riemann_compare.solvers_2D)
