#!/usr/bin/env python
# encoding: utf-8
r"""
Pyclaw testing framework

:Authors:
	Kyle Mandli (2010-09-03) Initial version
"""
# ============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import unittest

from test_evolve import *
from test_plotters import *
from test_util import *

__all__ = ['solver_test','util_test','plotter_test']

def run(test_suite='all',verbosity=2,debug=False):
    r"""Run the tests specified by test_suite
    
    :Input:
     - test_suite - (str) String representing the suite of tests to run, 
      available tests are,
       - 'solver' - Run all pyclaw solver tests
       - 'plotting' - Run all pyclaw plotting tests
       - 'util' - Run all utility function tests
       - 'all' - Run all test suites available, default option
    """
    # Create test suites
    # Solver tests
    solver_suite = unittest.TestSuite()
    solver_suite.addTest(LimiterTest())
    solver_suite.addTest(AdvectionTest())
    solver_suite.addTest(AcousticsTest())
    solver_suite.addTest(BurgersTest())
    solver_suite.addTest(EulerTest())
    solver_suite.addTest(ShallowTest())
    
    # Plotting tests
    plotter_suite = unittest.TestSuite()
    
    # Utility function tests
    data_manip_suite = unittest.TestLoader().loadTestsFromTestCase(DataManipulationTestCase)
    util_suite = unittest.TestSuite([data_manip_suite])

    if test_suite == 'solver':
        suite = solver_suite
    elif test_suite == 'plotting':
        suite = plotter_suite
    elif test_suite == 'util':
        suite = util_suite
    elif test_suite == 'all':
        suite = unittest.TestSuite([solver_suite,plotter_suite,util_suite])
    
    # Run suites
    if not debug:
        unittest.TextTestRunner(verbosity=verbosity).run(suite)
    else:
        suite.debug()
