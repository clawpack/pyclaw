import ctypes
import os.path

class iso_c_rp1_advection():

    def __init__(self, u):
        this_path = os.path.dirname(__file__)
        self._dll = ctypes.CDLL(os.path.join(this_path, '_iso_c_advection.so'))
        self.context = ctypes.c_void_p()
        self._u = ctypes.c_double(u)
        self._dll.rp1_advection_new(ctypes.byref(self._u),
                                   ctypes.byref(self.context))
        self.num_eqn = 1
        self.num_waves = 1
        self._rp = self._dll.rp1_advection_c

    def __del__(self):
        self._dll.rp1_advection_delete(ctypes.byref(self.context))
