'''
Created on Mar 18, 2012

@author: kristof
'''
import clawpack.pyclaw as pyclaw

class SubgridSolver(object):
    r"""
    The subgrid solver holds all information needed for the PyClaw/Clawpack solver
    to work on a single patch. It has to be thread safe to allow easy shared memory
    parallelization.
    
     
    """
    
    def __init__(self, solver, global_state, q, qbc, aux, position, size, subdivision_factor_x0, subdivision_factor_x1, subdivision_factor_x2, unknowns_per_cell, aux_fields_per_cell):
        r"""
        Initializes this subgrid solver. It get's all information to prepare a domain and state for a
        single subgrid. 
        
        :Input:
         -  *solver* - (:class:`pyclaw.Solver`) The PyClaw-solver used for advancing this subgrid in time.
         -  *global_state* - (:class:`pyclaw.State`) The global state. This is not the state used for
                            for the actual solving of the timestep on the subgrid.
         -  *q* - The array storing the current solution.
         -  *qbc* - The array storing the solution of the last timestep including the ghostlayer.
         -  *position* - A d-dimensional tuple holding the position of the grid in the computational domain.
                         This measures the continuous real-world position in floats, not the integer position 
                         in terms of cells. 
         -  *size* - A d-dimensional tuple holding the size of the grid in the computational domain. This
                     measures the continuous real-world size in floats, not the size in terms of cells.
         -  *subdivision_factor* - The number of cells in one dimension of this subgrid. At the moment only
                                    square subgrids are allowed, so the total number of cells in the subgrid
                                    (excluding the ghostlayer) is subdivision_factor x subdivision_factor.
         -  *unknowns_per_cell* - The number of equations or unknowns that are stored per cell of the subgrid.
        
        """
        self.solver = solver
        self.dim_x = pyclaw.Dimension('x',position[0],position[0] + size[0],subdivision_factor_x0)
        self.dim_y = pyclaw.Dimension('y',position[1],position[1] + size[1],subdivision_factor_x1)
        domain = pyclaw.Domain([self.dim_x,self.dim_y])
        subgrid_state = pyclaw.State(domain, unknowns_per_cell, aux_fields_per_cell)
        subgrid_state.q = q
        subgrid_state.aux = aux
        subgrid_state.problem_data = global_state.problem_data
        self.solution = pyclaw.Solution(subgrid_state, domain)

        
        self.solver.bc_lower[0] = pyclaw.BC.custom
        self.solver.bc_upper[0] = pyclaw.BC.custom
        self.solver.bc_lower[1] = pyclaw.BC.custom
        self.solver.bc_upper[1] = pyclaw.BC.custom
        
        self.qbc = qbc
        self.solver.user_bc_lower = self.user_bc_lower
        self.solver.user_bc_upper = self.user_bc_upper
        
    def step(self, maximum_timestep_size, estimated_next_dt):
        r"""
        Performs one timestep on the subgrid. This might result in several runs of the
        solver to find the maximum allowed timestep size in terms of stability.
        
        :Input:
         -  *maximum_timestep_size* - This is the maximum allowed timestep size in terms
                                      of the grid topology and the global timestep. I.e. 
                                      neighboring subgrids might forbid a timestep on this 
                                      subgrid. Also this subgrid is not allowed to advance 
                                      further than the the global timestep.
         -  *estimated_next_dt*- This is the estimation for the maximum allowed timestep size
                                 in terms of stability and results from the cfl number of the
                                 last timestep performed on this grid.
        """
        self.solver.dt = min(maximum_timestep_size, estimated_next_dt)
        self.solver.setup(self.solution)
        # Set qbc and timestep for the current patch
        self.solver.qbc = self.qbc
        self.solver.dt_max = maximum_timestep_size
        self.solver.evolve_to_time(self.solution)
        
        return self.solution.state.q
        
        
    def user_bc_lower(self, grid,dim,t,qbc,mbc):
        #Todo 3D: Rewrite for arbitrary dimensions
        if dim == self.dim_x:
            for i in range(mbc):
                qbc[:,:,i] = self.qbc[:,:,i]
        else:
            for i in range(mbc):
                qbc[:,i,:] = self.qbc[:,i,:]
        
    def user_bc_upper(self, grid,dim,t,qbc,mbc):
        #Todo 3D: Rewrite for arbitrary dimensions
        if dim == self.dim_x:
            for i in range(mbc):
                qbc[:,:,self.dim_y.num_cells+mbc+i] = self.qbc[:,:,self.dim_y.num_cells+mbc+i]
        else:
            for i in range(mbc):
                qbc[:,self.dim_x.num_cells+mbc+i,:] = self.qbc[:,self.dim_x.num_cells+mbc+i,:]
        
