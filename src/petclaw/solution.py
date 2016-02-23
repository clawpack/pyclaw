from clawpack import pyclaw
from clawpack.pyclaw.solution import Solution
import numpy as np

class Solution(Solution):
    """ Parallel Solution class.
    """
    __doc__ += pyclaw.util.add_parent_doc(pyclaw.Solution)

    def get_read_func(self, file_format):
        from clawpack.petclaw import io
        if file_format == 'petsc':
            return io.petsc.read
        elif file_format == 'hdf5':
            return io.hdf5.read
        else:
            raise ValueError("File format %s not supported." % file_format)

    def get_write_func(self, file_format):
        from clawpack.petclaw import io
        if 'petsc' in file_format:
            return io.petsc.write
        elif 'hdf5' in file_format:
            return io.hdf5.write
        else:
            raise ValueError("File format %s not supported." % file_format)

    def _init_ds_solution(self):
        """
        Initializes a downsampled version of the solution
        """
        import clawpack.petclaw as pyclaw

        ds_domain = pyclaw.Domain([pyclaw.Dimension(self.domain.patch.lower_global[i],
                                                    self.domain.patch.upper_global[i],
                                                    self.domain.patch.num_cells_global[i]/self.downsampling_factors[i],
                                                    self.domain.grid.dimensions[i].name)
                                       for i in range(self.domain.num_dim)], proc_sizes=self.state.q_da.getProcSizes())
        ds_state = self._init_ds_state(self.state)
        self._ds_solution = pyclaw.Solution(ds_state, ds_domain)
        self._ds_solution.t = self.t

    def _init_ds_state(self, state):
        """
        Returns a downsampled version of the given state object
        """
        import clawpack.petclaw as pyclaw

        ds_domain = pyclaw.Domain([pyclaw.Dimension(self.domain.patch.lower_global[i],
                                                    self.domain.patch.upper_global[i],
                                                    self.domain.patch.num_cells_global[i]/self.downsampling_factors[i],
                                                    self.domain.grid.dimensions[i].name)
                                   for i in range(self.domain.num_dim)], proc_sizes=state.q_da.getProcSizes())
        ds_state = pyclaw.State(ds_domain,state.num_eqn,state.num_aux)

        return ds_state

    def downsample(self, write_aux,write_p):
        """
        Returns a downsampled version of the solution by local averaging over the downsampling factors
        """
        for i  in range(len(self.states)):
            state = self.states[i]
            if i > 0:
                self.ds_solution.states.append(self._init_ds_state(self.downsample_factors, state))
            ds_state = self.ds_solution.states[i]

            # Downsample q
            if write_p:
                ds_state.p = self._downsample_global_to_local(state.get_q_da_for_ds(np.max(self.downsampling_factors)),
                                                             state.gqVec, state.num_eqn, ds_state.patch._da.getRanges())
            else:
                ds_state.q = self._downsample_global_to_local(state.get_q_da_for_ds(np.max(self.downsampling_factors)),
                                                             state.gqVec, state.num_eqn, ds_state.patch._da.getRanges())
            # Dowsample aux
            if write_aux:
                ds_state.aux = self._downsample_global_to_local(state.get_aux_da_for_ds(np.max(self.downsampling_factors)),
                                                               state.gauxVec, state.num_aux, ds_state.patch._da.getRanges())

        return self.ds_solution

    def _downsample_global_to_local(self, da_for_ds, gVec, num_eqn, ds_domain_ranges):
        """
        Returns a locally averaged array and handles ranges mapping between the original & downsampled solution objects

        Input:
            - da_for_ds: the DA object of the original solution (q or aux) with stencil width adjusted as appropriate to
              handle the boundaries of the array to be averaged because after domain decomposition, there are no guarantees that
              the downsampling factors will evenly divide the the sub-domain in each direction
            - gVec: The global vector associated with q or aux
            - num_eqn: number of components of q or aux
            - ds_domain_ranges: global ranges of the downsampled solution
        """

        from skimage.transform import downscale_local_mean
        """
        downscale_local_mean downsamples n-dimensional array by local averaging.
        First, it views the array as blocks of the downsampling factors, then it computes the local average of each block.

        Examples
        --------
        >>> a = np.arange(15).reshape(3, 5)
        >>> a
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14]])
        >>> downscale_local_mean(a, (2, 3))
        array([[ 3.5,  4. ],
               [ 5.5,  4.5]])

        """

        # Create local array with ghost cells
        q_local_vec = da_for_ds.createLocalVec()
        da_for_ds.globalToLocal(gVec, q_local_vec)
        q_local_array = q_local_vec.array
        shape = [end-start for start,end in da_for_ds.getGhostRanges()]
        shape.insert(0, num_eqn)
        q_local_array = q_local_array.reshape(shape, order='f')

        # Map global ranges in the original DA object to local ranges for local averaging
        q_da_ghost_ranges = da_for_ds.getGhostRanges()
        new_global_ranges = tuple([(s * self.downsampling_factors[i], e * self.downsampling_factors[i]) for i, (s,e) in enumerate(ds_domain_ranges)])
        new_local_slices = (slice(0, q_local_array.shape[0]),)
        new_local_slices = new_local_slices + tuple([slice(s - gs,q_local_array.shape[i+1] - (ge - e)) for i, ((s,e), (gs, ge)) in enumerate(zip(new_global_ranges, q_da_ghost_ranges))])
        q_local_array = q_local_array[new_local_slices]

        # Compute local mean
        return downscale_local_mean(q_local_array, (1,) + self.downsampling_factors)

