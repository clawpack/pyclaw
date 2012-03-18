'''
Created on Feb 7, 2012

@author: Kristof
'''
from pyclaw.solver import Solver
import signal
import logging
from ctypes import CDLL
from ctypes import c_double
from ctypes import c_int
from ctypes import c_void_p
from ctypes import CFUNCTYPE
from ctypes import py_object
from ctypes import POINTER

import pyclaw

class Solver(Solver):
    
    #Callback definitions
    CALLBACK_SOLVER = CFUNCTYPE(None, POINTER(c_double), py_object, py_object, c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double)
    CALLBACK_BOUNDARY_CONDITIONS = CFUNCTYPE(None, py_object, py_object, c_int, c_int)
    
    def __init__(self, solver, initialMinimalMeshWidth):
        self.solver = solver
        self.initialMinimalMeshWidth = initialMinimalMeshWidth
        self.dt_initial = solver.dt_initial
        self.num_ghost = solver.num_ghost
        self.rp = solver.rp
        
        self.solver_callback = self.get_solver_callback()
        self.boundary_condition_callback = self.get_boundary_condition_callback()
        
    def get_solver_callback(self):
        r"""
        Creates a closure for the solver callback method.
        """
        def callback_solver(return_dt_and_estimated_next_dt, q, qbc, subdivision_factor, unknowns_per_subcell, size, position_x, position_y, current_time, maximum_timestep_size, estimated_next_dt):
            # Set up grid information for current patch
            global dim_x
            dim_x = pyclaw.Dimension('x',position_x,position_x + size,subdivision_factor)
            global dim_y
            dim_y = pyclaw.Dimension('y',position_y,position_y + size,subdivision_factor)
            domain = pyclaw.Domain([dim_x,dim_y])
            state = pyclaw.State(domain, unknowns_per_subcell)
            state.problem_data = self.solution.state.problem_data
            state.q = q
            solution = pyclaw.Solution(state, domain)
            
            self.solver.bc_lower[0] = pyclaw.BC.custom
            self.solver.bc_upper[0] = pyclaw.BC.custom
            self.solver.bc_lower[1] = pyclaw.BC.custom
            self.solver.bc_upper[1] = pyclaw.BC.custom
            
            global ghostlayerArray 
            ghostlayerArray = qbc
            self.solver.user_bc_lower = self.user_bc_lower
            self.solver.user_bc_upper = self.user_bc_upper
            
            self.solver.dt = min(maximum_timestep_size, estimated_next_dt)
            
            self.solver.setup(solution)
            # Set qbc and timestep for the current patch
            self.solver.qbc = qbc
            self.solver.dt_max = maximum_timestep_size
            self.solver.evolve_to_time(solution)
            # Copy back the array with new values
            q[:]= state.q[:]
            
            return_dt_and_estimated_next_dt[0] = self.solver.dt
            if self.solver.cfl.get_cached_max() > 0:
                return_dt_and_estimated_next_dt[1] = self.solver.dt * self.solver.cfl_desired / self.solver.cfl.get_cached_max()
            else:
                return_dt_and_estimated_next_dt[1] = self.solver.dt_max
        return self.CALLBACK_SOLVER(callback_solver)
        
    def get_boundary_condition_callback(self):
        r"""
        Creates a closure for the boundary condition callback method.
        """
        def callback_boundary_conditions(q, qbc, dimension, setUpper):
            import numpy
            if(setUpper == 1):
                self.qbc_upper(self.solution.state, self.solution.state.grid.dimensions[dimension], self.solution.state.t, numpy.rollaxis(qbc,dimension+1,1), dimension)
            else:
                self.qbc_lower(self.solution.state, self.solution.state.grid.dimensions[dimension], self.solution.state.t, numpy.rollaxis(qbc,dimension+1,1), dimension)
        return self.CALLBACK_BOUNDARY_CONDITIONS(callback_boundary_conditions)
    
    def user_bc_lower(self, grid,dim,t,qbc,mbc):
#        print "Setting lower bc with mbc=" + str(mbc)
        if dim == dim_x:
#            print "Setting lower bc for x"
            for i in range(mbc):
    #            print "Setting " + str(boundaryArray[:,:,i])
                qbc[:,:,i] = ghostlayerArray[:,:,i]
        else:
#            print "Setting lower bc for y"
            for i in range(mbc):
    #            print "Setting " + str(boundaryArray[:,i,:])
                qbc[:,i,:] = ghostlayerArray[:,i,:]
    #    print "set qbc=" 
    #    print str(qbc)
        
    def user_bc_upper(self, grid,dim,t,qbc,mbc):
        if dim == dim_x:
#            print "Setting upper bc for x" + str(qbc.shape) + " " + str(ghostlayerArray.shape) + " mbc=" + str(mbc) + " dim=" + str(dim)
            for i in range(mbc):
    #            print "Setting (" +str(dim.n+mbc+i) + ") " + str(boundaryArray[:,:,dim.n+mbc+i])
                qbc[:,:,dim.num_cells+mbc+i] = ghostlayerArray[:,:,dim.num_cells+mbc+i]
        else:
#            print "Setting upper bc for y"
            for i in range(mbc):
    #            print "Setting (" +str(dim.n+mbc+i) + ") " + str(boundaryArray[:,dim.n+mbc+i,:])
                qbc[:,dim.num_cells+mbc+i,:] = ghostlayerArray[:,dim.num_cells+mbc+i,:]
    #    print "set qbc="
    #    print str(qbc)
    
    def setup(self, solution):
        logging.getLogger('peanoclaw').info("Loading Peano-library...")
        self.libpeano = CDLL('libpeano-claw-2d.dylib')
        logging.getLogger('peanoclaw').info("Peano loaded successfully.")
        self.libpeano.pyclaw_peano_new.restype = c_void_p
        self.libpeano.pyclaw_peano_destroy.argtypes = [c_void_p]
        self.libpeano.pyclaw_peano_evolveToTime.argtypes = [c_double, c_void_p, c_void_p, c_void_p]
        
        self.bc_lower = self.solver.bc_lower[:]
        self.bc_upper = self.solver.bc_upper[:]
        
        # Get parameters for Peano
        dimensions = solution.state.grid.dimensions
        subdivisionFactor = solution.state.grid.dimensions[0].num_cells
        
        self.libpeano.pyclaw_peano_new.argtypes = [c_double, c_double, c_double, c_int, c_double, c_void_p, c_void_p]
        self.peano = self.libpeano.pyclaw_peano_new(c_double(self.initialMinimalMeshWidth), \
                                                    c_double(dimensions[0].upper - dimensions[0].lower), \
                                                    c_double(dimensions[1].upper - dimensions[1].lower), \
                                                    subdivisionFactor,
                                                    self.solver.dt_initial,
                                                    self.boundary_condition_callback,
                                                    self.solver_callback)
        self.solver.setup(solution)
        
        # Set PeanoSolution
        import peanoclaw
        if(isinstance(solution, peanoclaw.Solution)):
            solution.peano = self.peano
            solution.libpeano = self.libpeano
        else:
            logging.getLogger('peanoclaw').warning("Use peanoclaw.Solution instead of pyclaw.Solution together with peanoclaw.Solver to provide plotting functionality.")
        
        #Causes Ctrl+C to quit Peano
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        
        self.libpeano.pyclaw_peano_runTests()
        
    def teardown(self):
        self.libpeano.pyclaw_peano_destroy(self.peano)
    
    def evolve_to_time(self,solution,tend=None):
        
        if(tend == None) :
            raise Exception("Not yet implemented.")
        
        self.solution = solution
        self.libpeano.pyclaw_peano_evolveToTime(tend, self.peano, self.boundary_condition_callback, self.solver_callback)
                
    def solve_one_timestep(self, q, qbc):
        self.solver.step(self.solution)
        