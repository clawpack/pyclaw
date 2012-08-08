'''
Created on Feb 7, 2012

@author: Kristof
'''
from clawpack.pyclaw.solver import Solver
import signal
import logging
from ctypes import CDLL
from ctypes import c_bool
from ctypes import c_double
from ctypes import c_int
from ctypes import c_void_p
from ctypes import c_char_p
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
    CALLBACK_INITIALIZATION = CFUNCTYPE(None, 
                                        py_object, #q
                                        py_object, #qbc
                                        py_object, #aux
                                        c_int,     #subdivision factor X0
                                        c_int,     #subdivision factor X1
                                        c_int,     #subdivision factor X2
                                        c_int,     #unknowns per cell
                                        c_int,     #aux fields per cell
                                        c_double, c_double, c_double, #size
                                        c_double, c_double, c_double) #position
    CALLBACK_SOLVER = CFUNCTYPE(c_double, 
                                POINTER(c_double), #Return array
                                py_object, #q
                                py_object, #qbc
                                py_object, #aux
                                c_int,     #subdivision factor X0
                                c_int,     #subdivision factor X1
                                c_int,     #subdivision factor X2
                                c_int,     #unknowns per cell
                                c_int,     #aux fields per cell
                                c_double, c_double, c_double, #size
                                c_double, c_double, c_double, #position
                                c_double, #current time
                                c_double, #maximum timestep size
                                c_double) #estimated next timestep size
    CALLBACK_BOUNDARY_CONDITIONS = CFUNCTYPE(None, py_object, py_object, c_int, c_int)
    
    def __init__(self, solver, initial_minimal_mesh_width, q_initialization, aux_initialization=None, refinement_criterion=None):
        r"""
        Initializes the Peano-solver. This keeps the Peano-spacetree internally and wraps the given PyClaw-solver.
        
        :Input:
         -  *solver* - (:class:`pyclaw.Solver`) The PyClaw-solver used internally.
         -  *initial_minimal_mesh_width* - The initial mesh width for the Peano mesh. I.e. Peano refines the mesh regularly
                                             until it is at least as fine as stated in this parameter.
        """
        self.solver = solver
        self.initial_minimal_mesh_width = initial_minimal_mesh_width
        self.q_initialization = q_initialization
        self.aux_initialization = aux_initialization
        self.refinement_criterion = refinement_criterion
        self.dt_initial = solver.dt_initial
        self.num_ghost = solver.num_ghost
        self.rp = solver.rp

        self.initialization_callback = self.get_initialization_callback()
        self.solver_callback = self.get_solver_callback()
        self.boundary_condition_callback = self.get_boundary_condition_callback()
        
    def get_initialization_callback(self):
        r"""
        Creates a closure for initializing the grid
        """
        def callback_initialization(q, qbc, aux, subdivision_factor_x0, subdivision_factor_x1, subdivision_factor_x2, unknowns_per_subcell, aux_fields_per_subcell, size_x, size_y, size_z, position_x, position_y, position_z):
            import clawpack.pyclaw as pyclaw
            self.dim_x = pyclaw.Dimension('x',position_x,position_x + size_x,subdivision_factor_x0)
            self.dim_y = pyclaw.Dimension('y',position_y,position_y + size_y,subdivision_factor_x1)
            #TODO 3D: use size_z and position_z
            domain = pyclaw.Domain([self.dim_x,self.dim_y])
            subgrid_state = pyclaw.State(domain, unknowns_per_subcell, aux_fields_per_subcell)
            subgrid_state.q = q
            if(aux_fields_per_subcell > 0):
              subgrid_state.aux = aux
            subgrid_state.problem_data = self.solution.state.problem_data
            self.q_initialization(subgrid_state)
            if(self.aux_initialization != None and aux_fields_per_subcell > 0):
              self.aux_initialization(subgrid_state)
            
        return self.CALLBACK_INITIALIZATION(callback_initialization)
        
    def get_solver_callback(self):
        r"""
        Creates a closure for the solver callback method.
        """
        def callback_solver(return_dt_and_estimated_next_dt, q, qbc, aux, subdivision_factor_x0, subdivision_factor_x1, subdivision_factor_x2, unknowns_per_cell, aux_fields_per_cell, size_x, size_y, size_z, position_x, position_y, position_z, current_time, maximum_timestep_size, estimated_next_dt):
        
            # Fix aux array
            if(aux_fields_per_cell == 0):
              aux = None
              
            #TODO 3D: Adjust position and size to 3D
            # Set up grid information for current patch
            import clawpack.peanoclaw as peanoclaw            
            subgridsolver = peanoclaw.SubgridSolver(self.solver, self.solution.state, q, qbc, aux, (position_x, position_y), (size_x, size_y), subdivision_factor_x0, subdivision_factor_x1, subdivision_factor_x2, unknowns_per_cell, aux_fields_per_cell)
            
            new_q = subgridsolver.step(maximum_timestep_size, estimated_next_dt)
            # Copy back the array with new values
            q[:]= new_q[:]
            
            return_dt_and_estimated_next_dt[0] = self.solver.dt
            if self.solver.cfl.get_cached_max() > 0:
              return_dt_and_estimated_next_dt[1] = self.solver.dt * self.solver.cfl_desired / self.solver.cfl.get_cached_max()
            else:
              return_dt_and_estimated_next_dt[1] = self.solver.dt_max
                
            #Steer refinement
            if not self.refinement_criterion == None:
            	return self.Trrefinement_criterion(subgridsolver.solution.state)
            else:
            	return self.initial_minimal_mesh_width
                
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
        self.libpeano = CDLL(self.get_lib_path())
        logging.getLogger('peanoclaw').info("Peano loaded successfully.")
        self.libpeano.pyclaw_peano_new.restype = c_void_p
        self.libpeano.pyclaw_peano_destroy.argtypes = [c_void_p]
        self.libpeano.pyclaw_peano_evolveToTime.argtypes = [c_double, c_void_p, c_void_p, c_void_p]
        
        self.bc_lower = self.solver.bc_lower[:]
        self.bc_upper = self.solver.bc_upper[:]
        self.user_bc_lower = self.solver.user_bc_lower
        self.user_bc_upper = self.solver.user_bc_upper
        
        # Get parameters for Peano
        dimensions = solution.state.grid.dimensions
        subdivision_factor_x0 = solution.state.grid.dimensions[0].num_cells
        subdivision_factor_x1 = solution.state.grid.dimensions[1].num_cells
        subdivision_factor_x2 = 0 #solution.state.grid.dimensions[2].num_cells #TODO 3D
        number_of_unknowns = solution.state.num_eqn 
        number_of_auxiliar_fields = solution.state.num_aux
        ghostlayer_width = self.num_ghost
        import os, sys
        configuration_file = os.path.join(sys.path[0], 'peanoclaw-config.xml')
        
        #self.solver.setup(solution)
        self.solution = solution
        
        self.libpeano.pyclaw_peano_new.argtypes = [ c_double, #Initial mesh width
                                                    c_double, #Domain position X0
                                                    c_double, #Domain position X1
                                                    c_double, #Domain position X2
                                                    c_double, #Domain size X0
                                                    c_double, #Domain size X1
                                                    c_double, #Domain size X2
                                                    c_int,    #Subdivision factor X0
                                                    c_int,    #Subdivision factor X1
                                                    c_int,    #Subdivision factor X2
                                                    c_int,    #Number of unknowns
                                                    c_int,    #Number of auxiliar fields
                                                    c_int,    #Ghostlayer width
                                                    c_double, #Initial timestep size
                                                    c_char_p, #Config file
                                                    c_bool,   #Run tests
                                                    c_void_p, #q Initialization callback
                                                    c_void_p, #Boundary condition callback
                                                    c_void_p] #Solver callback
        self.peano = self.libpeano.pyclaw_peano_new(c_double(self.initial_minimal_mesh_width),
                                                    c_double(dimensions[0].lower),
                                                    c_double(dimensions[1].lower),
                                                    c_double(0.0), #Todo 3D
                                                    c_double(dimensions[0].upper - dimensions[0].lower),
                                                    c_double(dimensions[1].upper - dimensions[1].lower),
                                                    c_double(0.0), #Todo 3D
                                                    subdivision_factor_x0,
                                                    subdivision_factor_x1,
                                                    subdivision_factor_x2,
                                                    number_of_unknowns,
                                                    number_of_auxiliar_fields,
                                                    ghostlayer_width,
                                                    self.solver.dt_initial,
                                                    c_char_p(configuration_file),
                                                    #True,
                                                    False,
                                                    self.initialization_callback,
                                                    self.boundary_condition_callback,
                                                    self.solver_callback)
        
        # Set PeanoSolution
        import clawpack.peanoclaw as peanoclaw
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
        self.solver.teardown()
    
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
        
    def get_lib_path(self):
        r"""
        Returns the path in which the shared library of Peano is located in.
        """
        import os
        import platform
        import clawpack.peanoclaw as peanoclaw
        if platform.system() == 'Linux':
            shared_library_extension = 'so'
        elif platform.system() == 'Darwin':
            shared_library_extension = 'dylib'
        else:
            raise("Unsupported operating system")    
        
        print(os.path.join(os.path.dirname(peanoclaw.__file__), 'libpeano-claw-2d.' + shared_library_extension))
        return os.path.join(os.path.dirname(peanoclaw.__file__), 'libpeano-claw-2d.' + shared_library_extension)
        

