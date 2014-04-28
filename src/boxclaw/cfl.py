
class CFL(object):
    """Parallel CFL object, responsible for computing the
    Courant-Friedrichs-Lewy condition across all processes.
    """

    def __init__(self, global_max):
        self._local_max = global_max
        self._global_max = global_max

    def get_global_max(self):
        r"""
        Compute the maximum CFL number over all processes for the current step.

        This is used to determine whether the CFL condition was
        violated and adjust the timestep.
        """
        import boxlib
        self._reduce_vec.array = self._local_max
        self._global_max = boxlib.bl[0].ReduceRealMax(self._local_max)
        return self._global_max

    def get_cached_max(self):
        return self._global_max

    def set_local_max(self,new_local_max):
        self._local_max = new_local_max

    def update_global_max(self,new_local_max):
        import boxlib
        self._global_max = boxlib.bl[0].ReduceRealMax(new_local_max)

