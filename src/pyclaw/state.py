r"""
Module containing all Pyclaw solution objects

:Authors:
    David I. Ketcheson -- Initial version (June 2011)
"""

import numpy as np

class State(object):
    r"""
    Contains the current state on a particular grid, including q, t, and aux.
    Usually aux is time-independent.

    Both q and aux are initialized to None.  They cannot be accessed until
    meqn and maux (respectively) are set.

    :State Data:
    
        The arrays :attr:`q`, and :attr:`aux` have variable 
        extents based on the grid dimensions and the values of 
        :attr:`meqn` and :attr:`maux`.  Note that these are initialy set to 
        None and later set to appropriately sized empty numpy arrays when
        :attr:`meqn` and :attr:`maux` are set.
 
    Typical usage, assuming a 1D grid::

        >>> state = State(grid)
        >>> meqn = 2
        >>> state.q[0,:] = sin(grid.x.center)
        >>> state.q[1,:] = cos(grid.x.center)
        >>> maux = 1
        >>> aux[0,:] = exp(grid.x.center)
    """

    # ========== Property Definitions ========================================
    def meqn():
        doc = r"""(int) - Number of unknowns (components of q)"""
        def fget(self):
            if self.q is None:
                raise Exception('state.meqn has not been set.')
            else: return self.q.shape[0]
        return locals()
    meqn = property(**meqn())

    def maux():
        doc = r"""(int) - Number of auxiliary fields"""
        def fget(self):
            if self.aux is not None: return self.aux.shape[0]
            else: return 0
        return locals()
    maux = property(**maux())

    def mp():
        doc = r"""(int) - Number of derived quantities"""
        def fset(self,mp):
            if self.p is not None:
                raise Exception('Cannot change state.mp after aux is initialized.')
            else:
                self.p = self.new_array(mp)
        def fget(self):
            if self.aux is not None: return self.aux.shape[0]
            else: return 0
        return locals()
    mp = property(**mp())

    def mF():
        doc = r"""(int) - Number of output functionals"""
        def fset(self,mF):
            if self.F is not None:
                raise Exception('Cannot change state.mF after aux is initialized.')
            else:
                self.F = self.new_array(mF)
        def fget(self):
            if self.F is not None: return self.F.shape[0]
            else: return 0
        return locals()
    mF = property(**mF())

    # ========== Class Methods ===============================================
    def __init__(self,grid,meqn,maux=0):
        import pyclaw.grid
        if not isinstance(grid,pyclaw.grid.Grid):
            raise Exception("""A PyClaw State object must be initialized with
                             a PyClaw Grid object.""")

        # ========== Attribute Definitions ===================================
        self.grid = grid
        r"""pyclaw.Grid.grid - The grid this state lives on"""
        self.p   = None
        r"""(ndarray(mp,...)) - Cell averages of derived quantities."""
        self.F   = None
        r"""(ndarray(mF,...)) - Cell averages of output functional densities."""
        self.aux_global = {}
        r"""(dict) - Dictionary of global values for this grid, 
            ``default = {}``"""
        self.t=0.
        r"""(float) - Current time represented on this grid, 
            ``default = 0.0``"""
        self.mcapa = -1

        self.q   = self.new_array(meqn)
        self.aux = self.new_array(maux)

    def __str__(self):
        output = "  t=%s meqn=%s\n  " % (self.t,self.meqn)
        if self.q is not None:
            output += "  q.shape=%s" % str(self.q.shape)
        if self.aux is not None:
            output += " aux.shape=%s" % str(self.aux.shape)
        return output

    def is_valid(self):
        r"""
        Checks to see if this state is valid
        
        The state is declared valid based on the following criteria:
            - :attr:`q` is not None
            - :attr:`meqn` > 0
            
        A debug logger message will be sent documenting exactly what was not 
        valid.
            
        :Output:
         - (bool) - True if valid, false otherwise.
        
        """
        import logging
        valid = True
        logger = logging.getLogger('solution')
        if self.q is None:
            logger.debug('The array q has not been initialized.')
            valid = False
        if self.meqn == 0:
            logger.debug('State.meqn has not been set.')
            valid = False
        if not self.q.flags['F_CONTIGUOUS']:
            logger.debug('q array is not Fortran contiguous.')
            valid = False
        return valid
 
    def set_cparam(self,fortran_module):
        """
        Set the variables in fortran_module.cparam to the corresponding values in
        grid.aux_global.  This is the mechanism for passing scalar variables to the
        Fortran Riemann solvers; cparam must be defined as a common block in the
        Riemann solver.

        This function should be called from solver.setup().  This seems like a fragile
        interdependency between solver and grid; perhaps aux_global should belong
        to solver instead of state.

        This function also checks that the set of variables defined in cparam 
        all appear in aux_global.
        """
        if hasattr(fortran_module,'cparam'):
            if not set(dir(fortran_module.cparam)) <= set(self.aux_global.keys()):
                raise Exception("""Some required value(s) in the cparam common 
                                   block in the Riemann solver have not been 
                                   set in aux_global.""")
            for global_var_name,global_var_value in self.aux_global.iteritems(): 
                setattr(fortran_module.cparam,global_var_name,global_var_value)

    def set_mbc(self,mbc):
        """
        Virtual routine (does nothing)
        """
        pass


    def set_q_from_qbc(self,mbc,qbc):
        """
        Set the value of q using the array qbc. for PetSolver, this
        involves setting qbc as the local vector array then perform
        a local to global communication. 
        """
        
        grid = self.grid
        if grid.ndim == 1:
            self.q = qbc[:,mbc:-mbc]
        elif grid.ndim == 2:
            self.q = qbc[:,mbc:-mbc,mbc:-mbc]
        elif grid.ndim == 3:
            self.q = qbc[:,mbc:-mbc,mbc:-mbc,mbc:-mbc]
        else:
            raise Exception("Assumption (1 <= ndim <= 3) violated.")

    def get_qbc_from_q(self,mbc,whichvec,qbc):
        """
        Fills in the interior of qbc (local vector) by copying q (global vector) to it.
        """
        ndim = self.grid.ndim
        
        if whichvec == 'q':
            q    = self.q
        elif whichvec == 'aux':
            q    = self.aux

        if ndim == 1:
            qbc[:,mbc:-mbc] = q
        elif ndim == 2:
            qbc[:,mbc:-mbc,mbc:-mbc] = q
        elif ndim == 3:
            qbc[:,mbc:-mbc,mbc:-mbc,mbc:-mbc] = q

        return qbc

    # ========== Copy functionality ==========================================
    def __copy__(self):
        return self.__class__(self)
        
        
    def __deepcopy__(self,memo={}):
        import copy
        result = self.__class__(copy.deepcopy(self.grid),self.meqn,self.maux)
        result.__init__(copy.deepcopy(self.grid),self.meqn,self.maux)
        
        for attr in ('t'):
            setattr(result,attr,copy.deepcopy(getattr(self,attr)))
        
        if self.q is not None:
            result.q = copy.deepcopy(self.q)
        if self.aux is not None:
            result.aux = copy.deepcopy(self.aux)
        result.aux_global = copy.deepcopy(self.aux_global)
        
        return result

    def sum_F(self,i):
        return np.sum(np.abs(self.F[i,...]))

    def new_array(self,dof):
        if dof==0: return None
        shape = [dof]
        shape.extend(self.grid.n)
        return np.empty(shape,order='F')
