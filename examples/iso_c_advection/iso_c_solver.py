from __future__ import absolute_import
import ctypes
import os.path
import numpy as np

from clawpack.pyclaw import ClawSolver1D

class iso_c_step1():

    def __init__(self):
        this_path = os.path.dirname(__file__)
        self._dll = ctypes.CDLL(os.path.join(this_path, 'iso_c_classic1.so'))

        # work arrays
        self.allocated = False
        self._f = None
        self._wave = None
        self._s = None
        self._amdq = None
        self._apdq = None
        self._dtdx = None

    def allocate(self, mx, num_eqn, num_ghost, num_waves):

        # these are work arrays, so order doesn't really matter
        self._f = np.empty((2*num_ghost + mx, num_eqn), dtype=np.double)
        self._wave = np.empty((num_eqn, num_waves, 2*num_ghost + mx),
                              dtype=np.double)
        self._s = np.empty((num_waves, 2*num_ghost+mx), dtype=np.double)
        self._amdq = np.empty((num_eqn, 2*num_ghost + mx), dtype=np.double)
        self._apdq = np.empty((num_eqn, 2*num_ghost + mx), dtype=np.double)
        self._dtdx = np.empty((2*num_ghost + mx), dtype=np.double)
        self.allocated = True

    def step1(self,
              py_num_ghost,
              py_mx,
              qbc,
              auxbc,
              py_dx,
              py_dt,
              method,
              mthlim,
              fwave,
              rp):
        r"""
        Take one time step on the homogeneous hyperbolic system.

        This function directly wraps the Clawpack step1 call, and is responsible
        for translating Pythonic data structures into their C/Fortran
        equivalents.
        """

        from ctypes import c_int, c_double, c_bool, byref

        # a real solver object would be caching/verifying these values, this is
        # just scaffolding

        py_num_eqn, mxbc = qbc.shape
        num_eqn = c_int(py_num_eqn)

        mthlim = np.asarray(mthlim, dtype=np.int)
        py_num_waves = mthlim.shape[0]
        num_waves = c_int(py_num_waves)

        py_num_aux, mxauxbc = auxbc.shape
        num_aux = c_int(py_num_aux)

        cfl = c_double()
        use_fwave = c_bool(False)

        num_ghost = c_int(py_num_ghost)
        mx = c_int(py_mx)

        dx = c_double(py_dx)
        dt = c_double(py_dt)

        if not self.allocated:
            self.allocate(py_mx, py_num_eqn, py_num_ghost, py_num_waves)


        def to_double_ref(nparray):
            return nparray.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


        def to_int_ref(nparray):
            return nparray.ctypes.data_as(ctypes.POINTER(ctypes.c_int))


        self._dll.step1_c(byref(num_eqn),
                          byref(num_waves),
                          byref(num_ghost),
                          byref(num_aux),
                          byref(mx),
                          to_double_ref(qbc),
                          to_double_ref(auxbc),
                          byref(dx),
                          byref(dt),
                          to_int_ref(method),
                          to_int_ref(mthlim),
                          byref(cfl),
                          to_double_ref(self._f),
                          to_double_ref(self._wave),
                          to_double_ref(self._s),
                          to_double_ref(self._amdq),
                          to_double_ref(self._apdq),
                          to_double_ref(self._dtdx),
                          byref(use_fwave),
                          byref(rp._rp),
                          byref(rp.context))


        return qbc, cfl.value


class ISO_C_ClawSolver1D(ClawSolver1D):

    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Create 1d ISO C Clawpack solver

        See :class:`ClawSolver1D` for more info.
        """

        self.iso_c_step1 = iso_c_step1()


        super(ISO_C_ClawSolver1D,self).__init__(riemann_solver,claw_package)


    def step_hyperbolic(self,solution):
        r"""
        Take one time step on the homogeneous hyperbolic system.

        :Input:
        - *solution* - (:class:`~pyclaw.solution.Solution`) Solution that
        will be evolved
        """

        state = solution.states[0]
        grid = state.grid

        self._apply_bcs(state)

        num_eqn,num_ghost = state.num_eqn,self.num_ghost

        mx = grid.num_cells[0]
        dx,dt = grid.delta[0],self.dt
        dtdx = np.zeros( (mx+2*num_ghost) ) + dt/dx

        self.qbc,cfl = self.iso_c_step1.step1(num_ghost, mx, self.qbc,
                                              self.auxbc, dx, dt, self._method,
                                              self._mthlim,self.fwave,
                                              self.rp)
        self.cfl.update_global_max(cfl)
        state.set_q_from_qbc(num_ghost,self.qbc)


        if state.num_aux > 0:
            state.set_aux_from_auxbc(num_ghost,self.auxbc)
