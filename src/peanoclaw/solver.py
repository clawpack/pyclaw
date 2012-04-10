'''
Created on Feb 7, 2012

@author: Kristof
'''
from pyclaw.solver import Solver
import signal
import logging
from ctypes import CDLL
from ctypes import c_bool
from ctypes import c_double
from ctypes import c_int
from ctypes import c_void_p
from ctypes import CFUNCTYPE
from ctypes import py_object
from ctypes import POINTER

class Solver(Solver):
    r"""
        This solver class wraps the AMR functionality of Peano. It holds a normal PyClaw-solver
        to advance separate subgrids in time. Therefore it provides the callbacks callback_solver(...)
        and callback_boundary_conditions(...) as an interface for Peano to use the PyClaw-solver.
        
        A Solver is typically instantiated as follows::

        >>> import pyclaw
        >>> solver = pyclaw.ClawSolver2D()
        >>> import peanoclaw
        >>> peanoclaw_solver = peanoclaw.Solver(solver, 1.0/18.0)
    """
    
    #Callback definitions
    CALLBACK_SOLVER = CFUNCTYPE(None, POINTER(c_double), py_object, py_object, c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double)
    CALLBACK_BOUNDARY_CONDITIONS = CFUNCTYPE(None, py_object, py_object, c_int, c_int)
    
    def __init__(self, solver, initial_minimal_mesh_width):
        r"""
        Initializes the Peano-solver. This keeps the Peano-spacetree internally and wraps the given PyClaw-solver.
        
        :Input:
         -  *solver* - (:class:`pyclaw.Solver`) The PyClaw-solver used internally.
         -  *initial_minimal_mesh_width* - The initial mesh width for the Peano mesh. I.e. Peano refines the mesh regularly
                                             until it is at least as fine as stated in this parameter.
        """
        self.solver = solver
        self.initial_minimal_mesh_width = initial_minimal_mesh_width
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
            import peanoclaw
            subgridsolver = peanoclaw.SubgridSolver(self.solver, self.solution.state, q, qbc, (position_x, position_y), (size, size), subdivision_factor, unknowns_per_subcell)
            
            new_q = subgridsolver.step(maximum_timestep_size, estimated_next_dt)
            # Copy back the array with new values
            q[:]= new_q[:]
            
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
    
    def setup(self, solution):
        r"""
        Initialize a Solver object. This method loads the library of Peano and prepares the initial mesh.
        
        See :class:`Solver` for full documentation
        """ 
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
        
        self.libpeano.pyclaw_peano_new.argtypes = [c_double, c_double, c_double, c_int, c_double, c_bool, c_void_p, c_void_p]
        self.peano = self.libpeano.pyclaw_peano_new(c_double(self.initial_minimal_mesh_width), \
                                                    c_double(dimensions[0].upper - dimensions[0].lower), \
                                                    c_double(dimensions[1].upper - dimensions[1].lower), \
                                                    subdivisionFactor,
                                                    self.solver.dt_initial,
                                                    True,
                                                    self.boundary_condition_callback,
                                                    self.solver_callback)
        self.solver.setup(solution)
        self.solution = solution
        
        # Set PeanoSolution
        import peanoclaw
        if(isinstance(solution, peanoclaw.Solution)):
            solution.peano = self.peano
            solution.libpeano = self.libpeano
        else:
            logging.getLogger('peanoclaw').warning("Use peanoclaw.Solution instead of pyclaw.Solution together with peanoclaw.Solver to provide plotting functionality.")
        
        #Causes Ctrl+C to quit Peano
        signal.signal(signal.SIGINT, signal.SIG_DFL)
                
    def teardown(self):
        r"""
        See :class:`Solver` for full documentation
        """ 
        self.libpeano.pyclaw_peano_destroy(self.peano)
    
    def evolve_to_time(self,solution,tend=None):
        r"""
        Performs one global timestep until all patches in the mesh reach the given end time.
        
        See :class:`Solver` for full documentation
        """ 
        if(tend == None) :
            raise Exception("Not yet implemented.")
        
        self.solution = solution
        self.libpeano.pyclaw_peano_evolveToTime(tend, self.peano, self.boundary_condition_callback, self.solver_callback)
                
    def solve_one_timestep(self, q, qbc):
        r"""
        """ 
        self.solver.step(self.solution)
        