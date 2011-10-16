#!/usr/bin/env python
# encoding: utf-8
r"""
Tests the evolve package

:Authors:
    Kyle Mandli (2010-09-04) Initial version
"""
# ============================================================================
#      Copyright (C) 2010 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import unittest
import os

import numpy as np

#import pyclaw.controller
#import pyclaw.evolve.clawpack
#import pyclaw.solution

_DATA_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),'data'))

# This is prevents nose from running these tests
from nose.plugins.skip import Skip,SkipTest
Skip()
raise SkipTest()

class SolverTestCase(unittest.TestCase):
    r"""Base solver test class"""
    
    def setUp(self):
        # Tolerance for comparison
        self.tolerance = 1e-6
        
        # Default settings
        self.controller = pyclaw.controller.Controller()
        self.controller.verbosity = 0
        self.controller.output_format = None
        self.controller.keep_copy = True
        self.controller.outstyle = 2
        self.controller.out_times = [0.0,1.0]
        self.controller.solution = self.true_solution(0.0)
        
        # Setup solver, using mainly defaults so should be setup later        
        self.controller.solver = pyclaw.evolve.clawpack.ClawSolver1D()
    
    def true_solution(self,t):
        r"""Returns the true or approximate solution at time t"""
        return None
        
    def runTest(self):
        # Run code
        self.controller.solution = self.true_solution(0.0)
        status = self.controller.run()
        
        # Compare each time output skipping the initial condition
        for sol in self.controller.frames[1:]:
            
            if self.controller.verbosity > 2:
                self.plot_solution()
            
            # Check the error for each frame
            self.assertTrue(np.all(self.error(sol,self.true_solution(sol.t)) 
                < self.tolerance))
            
    def error(self,solution,true_solution,index=None):
        r"""Error metric used when comparing solutions.
        
        """
        if index is None:
            index = range(solution.grids[0].q.shape[0])
        return np.linalg.norm(solution.grids[0].q[index,...] 
                    - true_solution.grids[0].q[index,...])
        
    def plot_solution(self,solution,markers=['k','xb','xg','xr','xc','xm'],
                            output=None,path="./",prefix="test",format="png"):
        r"""Plot the solution."""
        
        import matplotlib.pyplot as plt
        
        plt.clf()
        for m in solution.grids[0].q.shape[0]:
            plt.figure(m)
            plt.hold(True)
            plt.plot(solution.grids[0].dimensions[0].center,
                     solution.grids[0].q[m,:],'%s%s' % markers[m])
        if output == None:
            plt.show()
        else:
            for m in solution.grids[0].q.shape[0]:
                plt.figure(m)
                name = '%s%s.%s' % (prefix,str(m).zfill(4),format)
                plt.savefig(os.path.join(path,name))
        
    def write_solution(self,solution,**kargs):
        r"""Write out the solution to specify a new base test comparisons.
        
        """
        for (frame,solution) in enumerate(solutions):
            solution.write(frame+1,kargs)
            

class AdvectionTest(SolverTestCase):
    r"""Tests the 1d advection solver"""
    
    def setUp(self):
        # Call super set up routine
        super(AdvectionTest,self).setUp()
        
        # Customize solver setup
        self.controller.solver.dt = 0.1
        self.controller.solver.max_steps = 5000
        self.controller.solver.set_riemann_solver('advection')
        self.controller.solver.order = 2
        self.controller.solver.mthlim = 4
        self.controller.solver.mwaves = 1
    
    def true_solution(self,t):
    
        if t == 0.0:
            # Data parameters
            u = 1.0
            beta = 100.0
            gamma = 0.0
            x0 = 0.29999999999999999
            x1 = 0.7
            x2 = 0.9
    
            # Create empty solution
            x = pyclaw.solution.Dimension('x',0.0,1.0,100,bc_lower=2,bc_upper=2)
            grid = pyclaw.solution.Grid(x)
            grid.empty_q()
            grid.meqn = 1
            grid.aux_global['u'] = u
            grid.t = t
    
            # Gaussian
            qg = (np.exp(-beta * ((x.center-u*grid.t)-x0)**2) 
                            * np.cos(gamma * ((x.center-u*grid.t) - x0)))

            # Step Function
            qs = ((x.center-u*grid.t) > x1)*1.0 - ((x.center-u*grid.t) > x2)*1.0
    
            grid.q[0,:] = qg + qs
            
            sol = pyclaw.solution.Solution(grid)
        else:
            fname = os.path.join(_DATA_PATH,'advection_test')
            sol = pyclaw.solution.Solution(0,path=fname)
            
        return sol
        
        
class LimiterTest(AdvectionTest):
    r"""Tests for all the limiters implemented in the evolve.limiter package
    
    :note:
        Currently limiters 21 and 22 are not reproducible.  Limiter 19 is
        reproducible but suspect.  Should file a bug report on this!
    """
    
    def runTest(self):
        
        import pyclaw.evolve.limiters as limiters
        
        self.frames = []
        
        # Temporary fix, these limiters are not comparable to their base
        limiters = range(23)
        limiters.remove(19)
        limiters.remove(21)
        limiters.remove(22)
        for limiter in limiters:
            self.controller.solution = self.true_solution(0.0)
            self.controller.solver.mthlim = limiter
            status = self.controller.run()
            
            computed_solution = pyclaw.solution.Solution(limiter,
                            path=os.path.join(_DATA_PATH,'limiter_test'))
                                   
            print limiter
            self.assertTrue(abs(self.error(self.controller.solution,
                computed_solution)) < self.tolerance)


class AcousticsTest(SolverTestCase):
    
    def setUp(self):
        # Call super set up routine
        super(AcousticsTest,self).setUp()
        
        # Customize solver setup
        self.controller.solver.set_riemann_solver('acoustics')
        self.controller.solver.mthlim = [3,3]
        self.controller.solver.order = 2
        self.controller.solver.dt = 0.0001
        self.controller.solver.max_steps = 500
        self.controller.solver.mwaves = 2
    
    def true_solution(self,t):

        # Initialize Grid
        x = pyclaw.solution.Dimension('x',-1.0,1.0,100)
        grid = pyclaw.solution.Grid([x])
        grid.meqn = 2
        grid.t = t
        grid.mbc = 2
        grid.bc_lower = 2
        grid.bc_upper = 2
        
        # Problem parameters
        beta = 100.0
        gamma = 0.0
        x0 = 0.3
        grid.aux_global['rho'] = 1.0
        grid.aux_global['K'] =  4.0
        grid.aux_global['cc'] = np.sqrt(grid.aux_global['K'] 
            / grid.aux_global['rho'])
        grid.aux_global['zz'] = grid.aux_global['cc'] * grid.aux_global['rho']

        # Initial Condition
        grid.zeros_q()
        grid.q[0,:] = np.exp(-beta * (x.center - x0)**2) \
            * np.cos(gamma * (x.center - x0)**2)

        return pyclaw.solution.Solution(grid)

# This test does not work yet as teh true solution for time other than 0.0        
class BurgersTest(SolverTestCase):
    
    def setUp(self):
        # Call super setup
        super(BurgersTest,self).setUp()
        
        # Customize solver setup
        self.controller.solver.set_riemann_solver('burgers')
        self.controller.solver.mthlim = 3
        self.controller.solver.order = 2
        self.controller.solver.dt = 0.0001
        self.controller.solver.max_steps = 5000
        self.controller.solver.mwaves = 1
        self.controller.out_times = [0.0,0.5,1.0]
        
    def true_solution(self,t):
        # Setup problem parameters
        ic = 3
        beta = 100.0
        gamma = 0.0
        x0 = 0.3
        x1 = 0.6
        x2 = 0.8
        efix = False
    
        # Create empty grid
        x = pyclaw.solution.Dimension('x',0.0,1.0,100)
        grid = pyclaw.solution.Grid(x)
        grid.t = t
        grid.aux_global['efix'] = efix
        grid.empty_q()
    
        # Gaussian
        qg = np.exp(-beta * (x.center-x0)**2) * np.cos(gamma * (x.center - x0))

        # Step Function
        qs = (x.center > x1) * 1.0 - (x.center > x2) * 1.0
    
        if ic == 1:
            grid.q[0,:] = qg
        elif ic == 2:
            grid.q[0,:] = qs
        elif ic == 3:
            grid.q[0,:] = qg + qs
        
        return pyclaw.solution.Solution(grid)
        
        
# This test does not work yet as teh true solution for time other than 0.0
class EulerTest(SolverTestCase):
    
    def setUp(self):
        # Call super setup
        super(EulerTest,self).setUp()
        
        # Customize solver setup
        self.controller.solver.set_riemann_solver('euler_roe')
        self.controller.solver.mthlim = [4,4,4]
        self.controller.solver.order = 2
        self.controller.solver.dt = 0.1
        self.controller.solver.max_steps = 1000
        self.controller.solver.mwaves = 3
        
    def true_solution(self,t):
        # Setup problem parameters
        sloc = 0.5
        rhol = 3.0
        ul = 0.0
        pl = 3.0
        rhor = 1.0
        ur = 0.0
        pr = 1.0
        gamma = 1.4
        gamma1 = gamma - 1.0
        efix = False
    
        # Create empty grid
        x = pyclaw.solution.Dimension('x',0.0,1.0,100)
        grid = pyclaw.solution.Grid([x])
        grid.meqn = 3
        grid.aux_global['efix'] = efix
        grid.aux_global['gamma'] = gamma
        grid.aux_global['gamma1'] = gamma - 1.0
    
        El = pl/gamma1 + 0.5*rhol*ul**2
        Er = pr/gamma1 + 0.5*rhor*ur**2
    
        grid.empty_q()
        grid.q[0,:] = (x.center < sloc) * rhol + (x.center >= sloc) * rhor
        grid.q[1,:] = (x.center < sloc) * rhol*ul + (x.center >= sloc) * rhor*ur
        grid.q[2,:] = (x.center < sloc) * El + (x.center >= sloc) * Er

        return pyclaw.solution.Solution(grid)
        
        
# This test does not work yet as the true solution for time other than 0.0
class ShallowTest(SolverTestCase):
    
    def setUp(self):
        # Call super setup
        super(ShallowTest,self).setUp()
        
        # Customize solver setup
        self.controller.solver.dt = 0.1
        self.controller.solver.max_steps = 500
        self.controller.solver.set_riemann_solver('shallow_roe')
        self.controller.solver.order = 2
        self.controller.solver.mthlim = [3,3]
        self.controller.solver.mwaves = 2
        
    def true_solution(self,t):
        # Setup problem parameters
        g = 1.0
        sloc = 0.0
        hl = 3.0
        ul = 0.0
        hr = 1.0
        ur = 0.0
        efix = False
    
        # Create empty grid
        x = pyclaw.solution.Dimension('x',-5.0,5.0,500)
        grid = pyclaw.solution.Grid(x)
        grid.aux_global['g'] = g
        grid.aux_global['efix'] = efix
        grid.meqn = 2
        grid.t = t
        
        # Initial data
        grid.empty_q()
        grid.q[0,:] = hl * (grid.p_center[0] <= sloc) + hr * (grid.p_center[0] > sloc)
        grid.q[1,:] = hl*ul * (grid.p_center[0] <= sloc) + hr*ur * (grid.p_center[0] > sloc)

        return pyclaw.solution.Solution(grid)
        
        
if __name__ == '__main__':
    # Customize unittest due to empty SolverTestCase objectdef run(test_suite,verbosity=2,debug=False):
    suite = unittest.TestSuite()
    suite.addTest(LimiterTest())
    suite.addTest(AdvectionTest())
    # suite.addTest(AcousticsTest())
    # suite.addTest(BurgersTest())
    # suite.addTest(EulerTest())
    # suite.addTest(ShallowTest())

    unittest.TextTestRunner().run(suite)
