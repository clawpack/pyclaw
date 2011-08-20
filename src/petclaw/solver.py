r"""
This module contains the most abstract parallel solver class, PetSolver.
All parallel solvers inherit from this class.

:Authors:
    Amal Alghamdi
    David Ketcheson
"""
import pyclaw.solver

class CFL(object):
    def __init__(self, global_max):
        self._local_max = global_max
        self._global_max = global_max
        
    def get_global_max(self):
        r"""
        Compute the maximum CFL number over all processes for the current step.

        This is used to determine whether the CFL condition was
        violated and adjust the timestep.
        """
        from petsc4py import PETSc
        cflVec = PETSc.Vec().createWithArray([self._local_max])
        self._global_max = cflVec.max()[1]
        return self._global_max

    def get_cached_max(self):
        return self._global_max

    def set_local_max(self,new_local_max):
        self._local_max = new_local_max

    def update_global_max(self,new_local_max):
        from petsc4py import PETSc
        cflVec = PETSc.Vec().createWithArray([new_local_max])
        self._global_max = cflVec.max()[1]


# ============================================================================
#  Generic PetClaw solver class
# ============================================================================
class PetSolver(pyclaw.solver.Solver):
    r"""
    Generic PetClaw solver
    
    All PetClaw solvers inherit from this base class.
    See superclass pyclaw.solver.Solver for documentation of attributes.
    """

    def __init__(self,data=None):
        super(PetSolver,self).__init__(data)
        self.cfl = CFL(self._default_attr_values['cfl_desired'])
    
    # ========== Boundary Conditions ==================================

