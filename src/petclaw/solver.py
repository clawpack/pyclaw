r"""
This module contains the most abstract parallel solver class, PetSolver.
All parallel solvers inherit from this class.

:Authors:
    Amal Alghamdi
    David Ketcheson
"""
import pyclaw.solver

# ============================================================================
#  Generic PetClaw solver class
# ============================================================================
class PetSolver(pyclaw.solver.Solver):
    r"""
    Generic PetClaw solver
    
    All PetClaw solvers inherit from this base class.
    See superclass pyclaw.solver.Solver for documentation of attributes.
    """
    
    # ========== Boundary Conditions ==================================
    def copy_global_to_local(self,state,whichvec):
        """
        Returns q with ghost cells attached.  For PetSolver,
        this means returning the local vector.  
        """
        shape = [n + 2*self.mbc for n in state.grid.n]
        
        if whichvec == 'q':
            state.q_da.globalToLocal(state.gqVec, state.lqVec)
            shape.insert(0,state.meqn)
            self.qbc=state.lqVec.getArray().reshape(shape, order = 'F')

        elif whichvec == 'aux':
            state.aux_da.globalToLocal(state.gauxVec, state.lauxVec)
            shape.insert(0,state.maux)
            self.auxbc=state.lauxVec.getArray().reshape(shape, order = 'F')
    
    def set_global_q(self,state,ghosted_q):
        """
        Set the value of q using the array ghosted_q. for PetSolver, this
        involves setting ghosted_q as the local vector array then perform
        a local to global communication. 
        """
        
        #state.lqVec.placeArray(ghosted_q)
        #state.q_da.localToGlobal(state.lqVec,state.gqVec)
        #state.lqVec.resetArray() # This call is required because placeArray is
                                 # intended to be temporarly placement
        grid = state.grid
        if grid.ndim == 1:
            state.q = ghosted_q[:,self.mbc:-self.mbc]
        elif grid.ndim == 2:
            mbc, mx, my = self.mbc, grid.n[0],grid.n[1]
            state.q=ghosted_q[:,mbc:mx+mbc,mbc:my+mbc]
        else:
            raise NotImplementedError("The case of 3D is not handled in "\
            +"this function yet")

    def cfl():
        r"""
        Compute the maximum CFL number over all processes.

        This is used to determine whether the CFL condition was
        violated and adjust the timestep.
        """
        def fget(self):
            from petsc4py import PETSc
            if self.dt_variable:
                cflVec = PETSc.Vec().createWithArray([self._cfl])
                return cflVec.max()[1]
            else:
                return self._cfl
        def fset(self,cflnum):
            self._cfl = cflnum
        return locals()
    cfl = property(**cfl())

    def allocate_rk_stages(self,solution):
        r"""We could eliminate this function and just use
        the version in pyclaw.solver.Solver, if we were willing to
        check there whether the solver is a PetSolver.  But this
        would mean putting parallel-aware code in PyClaw, so for
        now we duplicate the function here.
        """
        from state import State

        if self.time_integrator   == 'Euler':  nregisters=1
        elif self.time_integrator == 'SSP33':  nregisters=2
        elif self.time_integrator == 'SSP104': nregisters=3
 
        state = solution.states[0]
        self._rk_stages = []
        for i in range(nregisters-1):
            self._rk_stages.append(State(state.grid,state.meqn,state.maux))
            self._rk_stages[-1].set_stencil_width(self.mbc)
            self._rk_stages[-1].aux_global       = state.aux_global
            self._rk_stages[-1].t                = state.t
            if state.maux > 0:
                self._rk_stages[-1].aux              = state.aux

