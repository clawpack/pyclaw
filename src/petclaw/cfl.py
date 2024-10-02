r"""
Module for the CFL object.
"""
class CFL(object):
    """ Parallel CFL object, responsible for computing the
    Courant-Friedrichs-Lewy condition across all processes.
    """

    def __init__(self, global_max):
        from petsc4py import PETSc
        self._local_max = global_max
        self._global_max = global_max
        self._reduce_vec = PETSc.Vec().createWithArray([0])

    def get_global_max(self):
        r"""
        Compute the maximum CFL number over all processes for the current step.
        This is used to determine whether the CFL condition was
        violated and adjust the timestep.
        """
        self._reduce_vec.array = self._local_max
        self._global_max = self._reduce_vec.max()[1]
        return self._global_max

    def get_cached_max(self):
        return self._global_max

    def set_local_max(self,new_local_max):
        self._local_max = new_local_max

    def set_global_max(self,new_global_max):
        self._global_max = new_global_max

    def update_global_max(self,new_local_max):
        self._reduce_vec.array = new_local_max
        self._global_max = max(self._global_max,self._reduce_vec.max()[1])
