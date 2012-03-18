'''
Created on Mar 17, 2012

@author: kristof
'''
import pyclaw
from ctypes import CFUNCTYPE
from ctypes import py_object
from ctypes import c_int
from ctypes import c_double
from ctypes import c_void_p

class Solution(pyclaw.solution.Solution):
    
    CALLBACK_ADD_PATCH_TO_SOLUTION = CFUNCTYPE(None, py_object, py_object, c_int, c_double, c_double, c_double, c_double)
    
    def __init__(self,*arg,**kargs):
        pyclaw.Solution.__init__(self,*arg,**kargs)
        self.peano = None
        self.libpeano = None
        
    def get_add_to_solution_callback(self):
        r"""
        Creates a closure for the callback method to add a grid to the solution.
        """
        def callback_add_to_solution(q, qbc, ghostlayer_width, size, position_X, position_Y, currentTime):
#            import pyclaw
            # Set up grid information for current patch
            subdivision_factor = q.shape[1]
            unknowns_per_subcell = q.shape[0]
            dim_x = pyclaw.Dimension('x', position_X, position_X + size, subdivision_factor)
            dim_y = pyclaw.Dimension('y', position_Y, position_Y + size, subdivision_factor)
#            domain = pyclaw.Domain([dim_x,dim_y])
#            state = pyclaw.State(domain, unknownsPerSubcell)
#            state.problem_data = self.solution.state.problem_data
#            state.q = q
#            solution = pyclaw.Solution(state, domain)

            patch = pyclaw.geometry.Patch((dim_x, dim_y))
            state = pyclaw.State(patch, unknowns_per_subcell)
            state.set_q_from_qbc(ghostlayer_width, qbc)
            
            print(patch)
            
            self.gathered_patches.append(patch)
            self.gathered_states.append(state)
            
        return self.CALLBACK_ADD_PATCH_TO_SOLUTION(callback_add_to_solution)

    
    def write(self,frame,path='./',file_format='ascii',file_prefix=None,
              write_aux=False,options={},write_p=False):
        if(self.peano == None or self.libpeano == None):
            raise Exception("self.peano and self.libpeano have to be initialized with the reference to the Peano grid by peanoclaw.Solver.setup(...)")

        #Gather patches and states        
        self.gathered_patches = []
        self.gathered_states = []
        
        self.libpeano.pyclaw_peano_gatherSolution.argtypes = [c_void_p, c_void_p]
        self.libpeano.pyclaw_peano_gatherSolution(self.peano, self.get_add_to_solution_callback())

        #Assemble solution and write file
        domain = pyclaw.Domain(self.gathered_patches)
        solution = pyclaw.Solution(self.gathered_states[1], domain)
        solution.write(frame, path, file_format,file_prefix,write_aux,options,write_p)
            
            
            
            