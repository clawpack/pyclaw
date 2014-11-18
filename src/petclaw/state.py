import clawpack.pyclaw

class FunctionSpace(clawpack.pyclaw.state.FunctionSpace):

    def __init__(self,patch,num_dof):
        super(FunctionSpace,self).__init__(patch,num_dof)
        self.da = self._create_DA()


    def _create_DA(self,num_ghost=0):
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

            DA = PETSc.DA().create(dim=self.patch.num_dim,
                                          dof=self.num_dof,
                                          sizes=self.patch.num_cells_global,
                                          periodic_type = periodic_type,
                                          stencil_width=num_ghost,
                                          comm=PETSc.COMM_WORLD)
        else:
            DA = PETSc.DA().create(dim=self.patch.num_dim,
                                          dof=self.num_dof,
                                          sizes=self.patch.num_cells_global,
                                          boundary_type = PETSc.DA.BoundaryType.PERIODIC,
                                          stencil_width=num_ghost,
                                          comm=PETSc.COMM_WORLD)

        return DA


class State(clawpack.pyclaw.State):
    """Parallel State class"""

    __doc__ += clawpack.pyclaw.util.add_parent_doc(clawpack.pyclaw.state)

    @property
    def num_eqn(self):
        r"""(int) - Number of unknowns (components of q)"""
        if self.q_space is None:
            raise Exception('state.num_eqn has not been set.')
        else: return self.q_space.num_dof

    @property
    def num_aux(self):
        r"""(int) - Number of auxiliary fields"""
        if self.aux_space is None: return 0
        else: return self.aux_space.num_dof

    @property
    def mp(self):
        r"""(int) - Number of derived quantities (components of p)"""
        if not hasattr(self,'p_space'):
            raise Exception('state.mp has not been set.')
        else: return self.p_space.num_dof
    @mp.setter
    def mp(self,mp):
        if hasattr(self,'p_space'):
            raise Exception('You cannot change state.mp after p is initialized.')
        else:
            self.p_space = FunctionSpace(self.patch,mp)
            self._p_global_vector = self.p_space.da.createGlobalVector()

    @property
    def mF(self):
        r"""(int) - Number of derived quantities (components of p)"""
        if not hasattr(self,'F_space'):
            raise Exception('state.mF has not been set.')
        else: return self.F_space.num_dof
    @mF.setter
    def mF(self,mF):
        if hasattr(self,'F_space'):
            raise Exception('You cannot change state.mF after F is initialized.')
        else:
            self.F_space = FunctionSpace(self.patch,mF)
            self._F_global_vector = self.F_space.da.createGlobalVector()

    @property
    def q(self):
        r"""
        Array of solution values.
        """
        shape = self.grid.num_cells
        shape.insert(0,self.num_eqn)
        return self._q_global_vector.getArray().reshape(shape, order = 'F')
    @q.setter
    def q(self,val):
        self._q_global_vector.setArray(val.reshape([-1], order = 'F'))

    @property
    def p(self):
        r"""
        Array containing values of derived quantities for output.
        """
        if not hasattr(self,'p_space'): return 0
        shape = self.grid.num_cells
        shape.insert(0,self.mp)
        p=self._p_global_vector.getArray().reshape(shape, order = 'F')
        return p
    @p.setter
    def p(self,val):
        mp = val.shape[0]
        if not hasattr(self,'_p_global_vector'): self.mp = mp
        self._p_global_vector.setArray(val.reshape([-1], order = 'F'))

    @property
    def F(self):
        r"""
        Array containing pointwise values (densities) of output functionals.
        This is just used as temporary workspace before summing.
        """
        if not hasattr(self,'F_space'): return 0
        shape = self.grid.num_cells
        shape.insert(0,self.mF)
        F=self._F_global_vector.getArray().reshape(shape, order = 'F')
        return F
    @F.setter
    def fset(self,val):
        mF = val.shape[0]
        if not hasattr(self,'_F_global_vector'): self.mF = mF
        self._F_global_vector.setArray(val.reshape([-1], order = 'F'))

    @property
    def aux(self):
        """
        We never communicate aux values; every processor should set its own ghost cell
        values for the aux array.  The global aux vector is used only for outputting
        the aux values to file; everywhere else we use the local vector.
        """
        if self.aux_space is None: return None
        shape = self.grid.num_cells
        shape.insert(0,self.num_aux)
        aux=self._aux_global_vector.getArray().reshape(shape, order = 'F')
        return aux
    @aux.setter
    def aux(self,val):
        # It would be nice to make this work also for parallel
        # loading from a file.
        if self.aux_space is None: 
            num_aux=val.shape[0]
            self.aux_space = FunctionSpace(num_aux)
        self._aux_global_vector.setArray(val.reshape([-1], order = 'F'))
    @property
    def num_dim(self):
        return self.patch.num_dim


    def __init__(self,geom,num_eqn,num_aux=0):
        r"""
        Here we don't call super because q and aux must be properties in PetClaw
        but should not be properties in PyClaw.

        The arguments num_eqn and num_aux can be integers, in which case
        a new DA with that many DOFs is created.  Alternatively, they
        can be functionSpaces, in which case the existing function space is
        just attached to the state.

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
                             a PetClaw Patch or Domain object.""")

        self.aux_space = None
        self.q_space = None

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

        if type(num_eqn) is int:
            self.q_space = FunctionSpace(self.patch,num_eqn)
        else:
            self.q_space = num_eqn

        if (type(num_aux) is int) and num_aux>0: 
            self.aux_space = FunctionSpace(self.patch,num_aux)
            self._aux_global_vector = self.aux_space.da.createGlobalVector()
        elif num_aux != 0:
            self.aux_space = num_aux
        else: # num_aux==0
            self.aux_space = None

        self._init_global_vecs()
        self._init_local_vecs()

    def _init_global_vecs(self):
        r"""
        Initializes PETSc global Vectors.
        """
        self._q_global_vector = self.q_space.da.createGlobalVector()
        if self.num_aux > 0:
            self._aux_global_vector = self.aux_space.da.createGlobalVector()
            
    def _init_local_vecs(self):
        r"""
        Initializes PETSc local Vectors.
        """
        self._q_local_vector = self.q_space.da.createLocalVector()
        if self.num_aux > 0:
            self._aux_local_vector = self.aux_space.da.createLocalVector()
 

    def get_qbc_from_q(self,num_ghost,qbc):
        """
        Returns q with ghost cells attached, by accessing the local vector.
        """
        shape = [n + 2*num_ghost for n in self.grid.num_cells]
        
        self.q_space.da.globalToLocal(self._q_global_vector, self._q_local_vector)
        shape.insert(0,self.num_eqn)
        return self._q_local_vector.getArray().reshape(shape, order = 'F')
            
    def get_auxbc_from_aux(self,num_ghost,auxbc):
        """
        Returns aux with ghost cells attached, by accessing the local vector.
        """
        shape = [n + 2*num_ghost for n in self.grid.num_cells]
        
        self.aux_space.da.globalToLocal(self._aux_global_vector, self._aux_local_vector)
        shape.insert(0,self.num_aux)
        return self._aux_local_vector.getArray().reshape(shape, order = 'F')

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
        self.q_space.da = self.q_space._create_DA(num_ghost)

        if self.num_aux > 0:
            aux0 = self.aux.copy()
            self.aux_space.da = self.aux_space._create_DA(num_ghost)

        # Need new vecs because we have new DAs
        self._init_global_vecs()
        self._init_local_vecs()

        # Copy old state into new vecs
        self.q = q0
        if self.num_aux > 0:
            self.aux = aux0

    def sum_F(self,i):
        return self._F_global_vector.strideNorm(i,0)

    def get_q_global(self):
        r"""
        Returns a copy of the global q array on process 0, otherwise returns None
        """
        from petsc4py import PETSc
        q_natural = self.q_space.da.createNaturalVec()
        self.q_space.da.globalToNatural(self._q_global_vector, q_natural)
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
        aux_natural = self.aux_space.da.createNaturalVec()
        self.aux_space.da.globalToNatural(self._aux_global_vector, aux_natural)
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
