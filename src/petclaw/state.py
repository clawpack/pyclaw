import clawpack.pyclaw

class State(clawpack.pyclaw.State):
    r"""  See the corresponding PyClaw class documentation."""

    @property
    def num_eqn(self):
        r"""(int) - Number of unknowns (components of q)"""
        if self.q_da is None:
            raise Exception('state.num_eqn has not been set.')
        else: return self.q_da.dof

    @property
    def mp(self):
        r"""(int) - Number of derived quantities (components of p)"""
        if self._p_da is None:
            raise Exception('state.mp has not been set.')
        else: return self._p_da.dof
    @mp.setter
    def mp(self,mp):
        if self._p_da is not None:
            raise Exception('You cannot change state.mp after p is initialized.')
        else:
            self._p_da = self._create_DA(mp)
            self.gpVec = self._p_da.createGlobalVector()

    @property
    def mF(self):
        r"""(int) - Number of derived quantities (components of p)"""
        if self._F_da is None:
            raise Exception('state.mF has not been set.')
        else: return self._F_da.dof
    @mF.setter
    def mF(self,mF):
        if self._F_da is not None:
            raise Exception('You cannot change state.mp after p is initialized.')
        else:
            self._F_da = self._create_DA(mF)
            self.gFVec = self._F_da.createGlobalVector()

    @property
    def num_aux(self):
        r"""(int) - Number of auxiliary fields"""
        if self.aux_da is None: return 0
        else: return self.aux_da.dof

    @property
    def q(self):
        r"""
        Array to store solution (q) values.

        Settting state.num_eqn automatically allocates space for q, as does
        setting q itself.
        """
        if self.q_da is None: return 0
        shape = self.grid.num_cells
        shape.insert(0,self.num_eqn)
        q=self.gqVec.getArray().reshape(shape, order = 'F')
        return q
    @q.setter
    def q(self,val):
        num_eqn = val.shape[0]
        if self.gqVec is None: self._init_q_da(num_eqn)
        self.gqVec.setArray(val.reshape([-1], order = 'F'))

    @property
    def p(self):
        r"""
        Array containing values of derived quantities for output.
        """
        if self._p_da is None: return 0
        shape = self.grid.num_cells
        shape.insert(0,self.mp)
        p=self.gpVec.getArray().reshape(shape, order = 'F')
        return p
    @p.setter
    def p(self,val):
        mp = val.shape[0]
        if self.gpVec is None: self.init_p_da(mp)
        self.gpVec.setArray(val.reshape([-1], order = 'F'))

    @property
    def F(self):
        r"""
        Array containing pointwise values (densities) of output functionals.
        This is just used as temporary workspace before summing.
        """
        if self._F_da is None: return 0
        shape = self.grid.num_cells
        shape.insert(0,self.mF)
        F=self.gFVec.getArray().reshape(shape, order = 'F')
        return F
    @F.setter
    def fset(self,val):
        mF = val.shape[0]
        if self.gFVec is None: self.init_F_da(mF)
        self.gFVec.setArray(val.reshape([-1], order = 'F'))

    @property
    def aux(self):
        """
        We never communicate aux values; every processor should set its own ghost cell
        values for the aux array.  The global aux vector is used only for outputting
        the aux values to file; everywhere else we use the local vector.
        """
        if self.aux_da is None: return None
        shape = self.grid.num_cells
        shape.insert(0,self.num_aux)
        aux=self.gauxVec.getArray().reshape(shape, order = 'F')
        return aux
    @aux.setter
    def aux(self,val):
        # It would be nice to make this work also for parallel
        # loading from a file.
        if self.aux_da is None: 
            num_aux=val.shape[0]
            self._init_aux_da(num_aux)
        self.gauxVec.setArray(val.reshape([-1], order = 'F'))
    @property
    def num_dim(self):
        return self.patch.num_dim


    def __init__(self,geom,num_eqn,num_aux=0):
        r"""
        Here we don't call super because q and aux must be properties in PetClaw
        but should not be properties in PyClaw.

        :attributes:
        patch - The patch this state lives on
        """

        from clawpack.pyclaw import geometry
        if isinstance(geom,geometry.Patch):
            self.patch = geom
        elif isinstance(geom,geometry.Domain):
            self.patch = geom.patches[0]
        else:
            raise Exception("""A PetClaw State object must be initialized with
                             a PyClaw Patch or Domain object.""")

        self.aux_da = None
        self.q_da = None

        self._p_da = None
        self.gpVec = None

        self._F_da = None
        self.gFVec = None

        # ========== Attribute Definitions ===================================
        self.problem_data = {}
        r"""(dict) - Dictionary of global values for this patch, 
            ``default = {}``"""
        self.t=0.
        r"""(float) - Current time represented on this patch, 
            ``default = 0.0``"""
        self.index_capa = -1
        self.keep_gauges = False
        r"""(bool) - Keep gauge values in memory for every time step, 
        ``default = False``"""
        self.gauge_data = []
        r"""(list) - List of numpy.ndarray objects. Each element of the list
        stores the values of the corresponding gauge if ``keep_gauges`` is set
        to ``True``"""

        self._init_q_da(num_eqn)
        if num_aux>0: self._init_aux_da(num_aux)

    def _init_aux_da(self,num_aux,num_ghost=0):
        r"""
        Initializes PETSc DA and global & local Vectors for handling the
        auxiliary array, aux. 
        
        Initializes aux_da, gauxVec and lauxVec.
        """
        self.aux_da = self._create_DA(num_aux,num_ghost)
        self.gauxVec = self.aux_da.createGlobalVector()
        self.lauxVec = self.aux_da.createLocalVector()
 
    def _init_q_da(self,num_eqn,num_ghost=0):
        r"""
        Initializes PETSc DA and Vecs for handling the solution, q. 
        
        Initializes q_da, gqVec and lqVec.
        """
        self.q_da = self._create_DA(num_eqn,num_ghost)
        self.gqVec = self.q_da.createGlobalVector()
        self.lqVec = self.q_da.createLocalVector()

    def _create_DA(self,dof,num_ghost=0):
        r"""Returns a PETSc DA and associated global Vec.
        Note that no local vector is returned.
        """
        from petsc4py import PETSc

        #Due to the way PETSc works, we just make the patch always periodic,
        #regardless of the boundary conditions actually selected.
        #This works because in solver.qbc() we first call globalToLocal()
        #and then impose the real boundary conditions (if non-periodic).

        if hasattr(PETSc.DA, 'PeriodicType'):
            if self.num_dim == 1:
                periodic_type = PETSc.DA.PeriodicType.X
            elif self.num_dim == 2:
                periodic_type = PETSc.DA.PeriodicType.XY
            elif self.num_dim == 3:
                periodic_type = PETSc.DA.PeriodicType.XYZ
            else:
                raise Exception("Invalid number of dimensions")

            DA = PETSc.DA().create(dim=self.num_dim,
                                          dof=dof,
                                          sizes=self.patch.num_cells_global,
                                          periodic_type = periodic_type,
                                          stencil_width=num_ghost,
                                          comm=PETSc.COMM_WORLD)
        else:
            DA = PETSc.DA().create(dim=self.num_dim,
                                          dof=dof,
                                          sizes=self.patch.num_cells_global,
                                          boundary_type = PETSc.DA.BoundaryType.PERIODIC,
                                          stencil_width=num_ghost,
                                          comm=PETSc.COMM_WORLD)

        return DA

    def set_q_from_qbc(self,num_ghost,qbc):
        """
        Set the value of q using the array qbc. for PetSolver, this
        involves setting qbc as the local vector array then perform
        a local to global communication. 
        """
        
        patch = self.patch
        if patch.num_dim == 1:
            self.q = qbc[:,num_ghost:-num_ghost]
        elif patch.num_dim == 2:
            self.q = qbc[:,num_ghost:-num_ghost,num_ghost:-num_ghost]
        elif patch.num_dim == 3:
            self.q = qbc[:,num_ghost:-num_ghost,num_ghost:-num_ghost,num_ghost:-num_ghost]
        else:
            raise NotImplementedError("The case of 3D is not handled in "\
            +"this function yet")

    def set_aux_from_auxbc(self,num_ghost,auxbc):
        """
        Set the value of aux using the array auxbc. for PetSolver, this
        involves setting auxbc as the local vector array then perform
        a local to global communication. 
        """
        
        patch = self.patch
        if patch.num_dim == 1:
            self.aux = auxbc[:,num_ghost:-num_ghost]
        elif patch.num_dim == 2:
            self.aux = auxbc[:,num_ghost:-num_ghost,num_ghost:-num_ghost]
        elif patch.num_dim == 3:
            self.aux = auxbc[:,num_ghost:-num_ghost,num_ghost:-num_ghost,num_ghost:-num_ghost]
        else:
            raise NotImplementedError("The case of 3D is not handled in "\
            +"this function yet")


    def get_qbc_from_q(self,num_ghost,qbc):
        """
        Returns q with ghost cells attached.  For PetSolver,
        this means returning the local vector.  
        """
        shape = [n + 2*num_ghost for n in self.grid.num_cells]
        
        self.q_da.globalToLocal(self.gqVec, self.lqVec)
        shape.insert(0,self.num_eqn)
        return self.lqVec.getArray().reshape(shape, order = 'F')
            
    def get_auxbc_from_aux(self,num_ghost,auxbc):
        """
        Returns aux with ghost cells attached.  For PetSolver,
        this means returning the local vector.  
        """
        shape = [n + 2*num_ghost for n in self.grid.num_cells]
        
        self.aux_da.globalToLocal(self.gauxVec, self.lauxVec)
        shape.insert(0,self.num_aux)
        return self.lauxVec.getArray().reshape(shape, order = 'F')

    def set_num_ghost(self,num_ghost):
        r"""
        This is a hack to deal with the fact that petsc4py
        doesn't allow us to change the stencil_width (num_ghost).

        Instead, we initially create DAs with stencil_width=0.
        Then, in solver.setup(), we call this function to replace
        those DAs with new ones that have the right stencil width.

        This could be made more efficient using some PETSc calls,
        but it only happens once so it seems not to be worth it.
        """
        q0 = self.q.copy()
        self._init_q_da(self.num_eqn,num_ghost)
        self.q = q0

        if self.aux is not None:
            aux0 = self.aux.copy()
            self._init_aux_da(self.num_aux,num_ghost)
            self.aux = aux0

    def sum_F(self,i):
        return self.gFVec.strideNorm(i,0)

    def get_q_global(self):
        r"""
        Returns a copy of the global q array on process 0, otherwise returns None
        """
        from petsc4py import PETSc
        q_natural = self.q_da.createNaturalVec()
        self.q_da.globalToNatural(self.gqVec, q_natural)
        scatter, q0Vec = PETSc.Scatter.toZero(q_natural)
        scatter.scatter(q_natural, q0Vec, False, PETSc.Scatter.Mode.FORWARD)
        rank = PETSc.COMM_WORLD.getRank()
        if rank == 0:
            shape = self.patch.num_cells_global
            shape.insert(0,self.num_eqn)
            q0=q0Vec.getArray().reshape(shape, order = 'F').copy()
        else:
            q0=None
        
        scatter.destroy()
        q0Vec.destroy()

        return q0

    def get_aux_global(self):
        r"""
        Returns a copy of the global aux array on process 0, otherwise returns None
        """
        from petsc4py import PETSc
        aux_natural = self.aux_da.createNaturalVec()
        self.aux_da.globalToNatural(self.gauxVec, aux_natural)
        scatter, aux0Vec = PETSc.Scatter.toZero(aux_natural)
        scatter.scatter(aux_natural, aux0Vec, False, PETSc.Scatter.Mode.FORWARD)
        rank = PETSc.COMM_WORLD.getRank()
        if rank == 0:
            shape = self.patch.num_cells_global
            shape.insert(0,self.num_aux)
            aux0=aux0Vec.getArray().reshape(shape, order = 'F').copy()
        else:
            aux0=None
        
        scatter.destroy()
        aux0Vec.destroy()

        return aux0
