
import fboxlib
import clawpack.pyclaw
import numpy as np

class State(clawpack.pyclaw.State):
    """Parallel State class"""

    def get_qbc_from_q(self,num_ghost,qbc):
        """
        Fills in the interior of qbc by copying q to it.
        """
        self._fill_boundary(self.num_eqn, num_ghost, self.q, qbc)
        return qbc

    def get_auxbc_from_aux(self,num_ghost,auxbc):
        self._fill_boundary(self.num_aux, num_ghost, self.aux, auxbc)
        return auxbc

    def _fill_boundary(self, ncomp, nghost, q, qbc):
        mf = self._create_multifab(ncomp, nghost)
        self._copy_into_multifab(mf, nghost, q)
        mf.fill_boundary()
        self._copy_outof_multifab(mf, qbc)
        del mf

    def _create_multifab(self, ncomp, nghost):
        return fboxlib.multifab(self.patch._la, ncomp, nghost)

    def _get_array(self, mf):
        return mf.fab(1).array

    def _copy_into_multifab(self, mf, ng, q):

        num_dim = self.patch.num_dim
        fab     = self._get_array(mf)

        if ng == 0:
            fab[...] = np.rollaxis(q, 0, q.ndim)
        elif num_dim == 1:
            fab[ng:-ng,:] = np.rollaxis(q, 0, q.ndim)
        elif num_dim == 2:
            fab[ng:-ng,ng:-ng,:] = np.rollaxis(q, 0, q.ndim)
        elif num_dim == 3:
            fab[ng:-ng,ng:-ng,ng:-ng,:] = np.rollaxis(q, 0, q.ndim)
        else:
            raise Exception("Assumption (1 <= num_dim <= 3) violated.")

    def _copy_outof_multifab(self, mf, q):
        fab = self._get_array(mf)
        q[...] = np.rollaxis(fab, self.q.ndim-1)
