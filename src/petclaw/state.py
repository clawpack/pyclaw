import pyclaw.state

class State(pyclaw.state.State):
    r"""  See the corresponding PyClaw class documentation."""

    def meqn():
        doc = r"""(int) - Number of unknowns (components of q)"""
        def fget(self):
            if self.q_da is None:
                raise Exception('state.meqn has not been set.')
            else: return self.q_da.dof
        return locals()
    meqn = property(**meqn())

    def mp():
        doc = r"""(int) - Number of derived quantities (components of p)"""
        def fset(self,mp):
            if self._p_da is not None:
                raise Exception('You cannot change state.mp after p is initialized.')
            else:
                self._p_da = self._create_DA(mp)
                self.gpVec = self._p_da.createGlobalVector()
        def fget(self):
            if self._p_da is None:
                raise Exception('state.mp has not been set.')
            else: return self._p_da.dof
        return locals()
    mp = property(**mp())

    def mF():
        doc = r"""(int) - Number of derived quantities (components of p)"""
        def fset(self,mF):
            if self._F_da is not None:
                raise Exception('You cannot change state.mp after p is initialized.')
            else:
                self._F_da = self._create_DA(mF)
                self.gFVec = self._F_da.createGlobalVector()
        def fget(self):
            if self._F_da is None:
                raise Exception('state.mF has not been set.')
            else: return self._F_da.dof
        return locals()
    mF = property(**mF())

    def maux():
        doc = r"""(int) - Number of auxiliary fields"""
        def fget(self):
            if self.aux_da is None: return 0
            else: return self.aux_da.dof
        return locals()
    maux = property(**maux())

    def q():
        r"""
        Array to store solution (q) values.

        Settting state.meqn automatically allocates space for q, as does
        setting q itself.
        """
        def fget(self):
            if self.q_da is None: return 0
            shape = self.grid.ng
            shape.insert(0,self.meqn)
            q=self.gqVec.getArray().reshape(shape, order = 'F')
            return q
        def fset(self,q):
            meqn = q.shape[0]
            if self.gqVec is None: self._init_q_da(meqn)
            self.gqVec.setArray(q.reshape([-1], order = 'F'))
        return locals()
    q = property(**q())

    def p():
        r"""
        Array containing values of derived quantities for output.
        """
        def fget(self):
            if self._p_da is None: return 0
            shape = self.grid.ng
            shape.insert(0,self.mp)
            p=self.gpVec.getArray().reshape(shape, order = 'F')
            return p
        def fset(self,p):
            mp = p.shape[0]
            if self.gpVec is None: self.init_p_da(mp)
            self.gpVec.setArray(p.reshape([-1], order = 'F'))
        return locals()
    p = property(**p())

    def F():
        r"""
        Array containing pointwise values (densities) of output functionals.
        This is just used as temporary workspace before summing.
        """
        def fget(self):
            if self._F_da is None: return 0
            shape = self.grid.ng
            shape.insert(0,self.mF)
            F=self.gFVec.getArray().reshape(shape, order = 'F')
            return F
        def fset(self,F):
            mF = F.shape[0]
            if self.gFVec is None: self.init_F_da(mF)
            self.gFVec.setArray(F.reshape([-1], order = 'F'))
        return locals()
    F = property(**F())

    def aux():
        """
        We never communicate aux values; every processor should set its own ghost cell
        values for the aux array.  The global aux vector is used only for outputting
        the aux values to file; everywhere else we use the local vector.
        """
        def fget(self):
            if self.aux_da is None: return None
            shape = self.grid.ng
            shape.insert(0,self.maux)
            aux=self.gauxVec.getArray().reshape(shape, order = 'F')
            return aux
        def fset(self,aux):
            # It would be nice to make this work also for parallel
            # loading from a file.
            if self.aux_da is None: 
                maux=aux.shape[0]
                self._init_aux_da(maux)
            self.gauxVec.setArray(aux.reshape([-1], order = 'F'))
        return locals()
    aux = property(**aux())
    def ndim():
        def fget(self):
            return self.grid.ndim
        return locals()
    ndim = property(**ndim())


    def __init__(self,grid,meqn,maux=0):
        r"""
        Here we don't call super because q and aux must be properties in PetClaw
        but should not be properties in PyClaw.

        """
        import petclaw.grid
        if not isinstance(grid,petclaw.grid.Grid):
            raise Exception("""A PetClaw State object must be initialized with
                             a PetClaw Grid object.""")
        self.aux_da = None
        self.q_da = None

        self._p_da = None
        self.gpVec = None

        self._F_da = None
        self.gFVec = None

        # ========== Attribute Definitions ===================================
        self.grid = grid
        r"""pyclaw.Grid.grid - The grid this state lives on"""
        self.aux_global = {}
        r"""(dict) - Dictionary of global values for this grid, 
            ``default = {}``"""
        self.t=0.
        r"""(float) - Current time represented on this grid, 
            ``default = 0.0``"""
        self.mcapa = -1

        self._init_q_da(meqn)
        if maux>0: self._init_aux_da(maux)

    def _init_aux_da(self,maux,mbc=0):
        r"""
        Initializes PETSc DA and global & local Vectors for handling the
        auxiliary array, aux. 
        
        Initializes aux_da, gauxVec and lauxVec.
        """
        self.aux_da = self._create_DA(maux,mbc)
        self.gauxVec = self.aux_da.createGlobalVector()
        self.lauxVec = self.aux_da.createLocalVector()
 
    def _init_q_da(self,meqn,mbc=0):
        r"""
        Initializes PETSc DA and Vecs for handling the solution, q. 
        
        Initializes q_da, gqVec and lqVec,
        and also sets up nstart, nend, and mbc for the dimensions.
        """
        self.q_da = self._create_DA(meqn,mbc)
        self.gqVec = self.q_da.createGlobalVector()
        self.lqVec = self.q_da.createLocalVector()

        #Now set the local indices for the Dimension objects:
        ranges = self.q_da.getRanges()
        for i,nrange in enumerate(ranges):
            dim = self.grid.dimensions[i]
            dim.nstart = nrange[0]
            dim.nend   = nrange[1]
            dim.lowerg = dim.lower + dim.nstart*dim.d

    def _create_DA(self,dof,mbc=0):
        r"""Returns a PETSc DA and associated global Vec.
        Note that no local vector is returned.
        """
        from petsc4py import PETSc

        #Due to the way PETSc works, we just make the grid always periodic,
        #regardless of the boundary conditions actually selected.
        #This works because in solver.qbc() we first call globalToLocal()
        #and then impose the real boundary conditions (if non-periodic).

        if hasattr(PETSc.DA, 'PeriodicType'):
            if self.ndim == 1:
                periodic_type = PETSc.DA.PeriodicType.X
            elif self.ndim == 2:
                periodic_type = PETSc.DA.PeriodicType.XY
            elif self.ndim == 3:
                periodic_type = PETSc.DA.PeriodicType.XYZ
            else:
                raise Exception("Invalid number of dimensions")

            DA = PETSc.DA().create(dim=self.ndim,
                                          dof=dof,
                                          sizes=self.grid.n,
                                          periodic_type = periodic_type,
                                          stencil_width=mbc,
                                          comm=PETSc.COMM_WORLD)
        else:
            DA = PETSc.DA().create(dim=self.ndim,
                                          dof=dof,
                                          sizes=self.grid.n,
                                          boundary_type = PETSc.DA.BoundaryType.PERIODIC,
                                          stencil_width=mbc,
                                          comm=PETSc.COMM_WORLD)

        return DA

    def set_q_from_qbc(self,mbc,qbc):
        """
        Set the value of q using the array qbc. for PetSolver, this
        involves setting qbc as the local vector array then perform
        a local to global communication. 
        """
        
        #state.lqVec.placeArray(ghosted_q)
        #self.q_da.localToGlobal(state.lqVec,state.gqVec)
        #state.lqVec.resetArray() # This call is required because placeArray is
                                 # intended to be temporarly placement
        grid = self.grid
        if grid.ndim == 1:
            self.q = qbc[:,mbc:-mbc]
        elif grid.ndim == 2:
            self.q = qbc[:,mbc:-mbc,mbc:-mbc]
        else:
            raise NotImplementedError("The case of 3D is not handled in "\
            +"this function yet")

    def get_qbc_from_q(self,mbc,whichvec,qbc):
        """
        Returns q with ghost cells attached.  For PetSolver,
        this means returning the local vector.  
        """
        shape = [n + 2*mbc for n in self.grid.ng]
        
        if whichvec == 'q':
            self.q_da.globalToLocal(self.gqVec, self.lqVec)
            shape.insert(0,self.meqn)
            return self.lqVec.getArray().reshape(shape, order = 'F')
            
        elif whichvec == 'aux':
            self.aux_da.globalToLocal(self.gauxVec, self.lauxVec)
            shape.insert(0,self.maux)
            return self.lauxVec.getArray().reshape(shape, order = 'F')

    def set_mbc(self,mbc):
        r"""
        This is a hack to deal with the fact that petsc4py
        doesn't allow us to change the stencil_width (mbc).

        Instead, we initially create DAs with stencil_width=0.
        Then, in solver.setup(), we call this function to replace
        those DAs with new ones that have the right stencil width.

        This could be made more efficient using some PETSc calls,
        but it only happens once so it seems not to be worth it.
        """
        q0 = self.q.copy()
        self._init_q_da(self.meqn,mbc)
        self.q = q0

        if self.aux is not None:
            aux0 = self.aux.copy()
            self._init_aux_da(self.maux,mbc)
            self.aux = aux0

    def sum_F(self,i):
        return self.gFVec.strideNorm(i,0)
