
import clawpack.pyclaw
import numpy as np

class State(clawpack.pyclaw.State):
    """Parallel State class"""

    def get_qbc_from_q(self,num_ghost,qbc):
        """
        XXX: Fills in the interior of qbc by copying q to it.
        """

        num_dim = self.patch.num_dim

        mf = self._create_multifab(self.num_eqn, num_ghost)
        a  = self._get_array(mf)

        if num_dim == 1:
            a[num_ghost:-num_ghost,:] = np.rollaxis(self.q, self.q.ndim-1)
        elif num_dim == 2:
            a[num_ghost:-num_ghost,num_ghost:-num_ghost,:] = np.rollaxis(self.q, self.q.ndim-1)
        elif num_dim == 3:
            a[num_ghost:-num_ghost,num_ghost:-num_ghost,num_ghost:-num_ghost,:] = np.rollaxis(self.q, self.q.ndim-1)
        else:
            raise Exception("Assumption (1 <= num_dim <= 3) violated.")

        mf.FillBoundary()
        self.patch._geom.FillPeriodicBoundary(mf)

        qbc[...] = np.rollaxis(a, self.q.ndim-1)

        return qbc


    def _create_multifab(self, ncomp, nghost):
        import boxlib
        return boxlib.MultiFab(self.patch._ba, ncomp, nghost)


    def _get_array(self, mf):
        for i in range(mf.size()):
            if mf[i] is not None:
                return mf[i].get_array()
        return None
