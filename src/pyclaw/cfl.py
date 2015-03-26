r"""
Module for the CFL object, which is responsible for computing and enforcing the
Courant-Friedrichs-Lewy condition.
"""

class CFL(object):
    def __init__(self, global_max):
        self._global_max = global_max
        
    def get_global_max(self):
        r"""
        Compute the maximum CFL number over all processes for the current step.

        This is used to determine whether the CFL condition was
        violated and adjust the timestep.
        """
        return self._global_max

    def get_cached_max(self):
        return self._global_max

    def set_local_max(self,new_local_max):
        self._global_max = new_local_max

    def set_global_max(self,new_global_max):
        self._global_max = new_global_max

    def update_global_max(self,new_local_max):
        self._global_max = max(self._global_max,new_local_max)
