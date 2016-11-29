r"""(Py)WENO based reconstructor for hyperbolic PDEs.

The :py:mod:`weno.reconstruct` module needs to be built before this
module can be used.  See 'weno/codegen.py' for details.

To build a higher order reconstruction, *k* needs to be tweaked here
and in 'weno/codegen.py'.  Also, *num_ghost* needs to be tweaked in the
PyClaw solver.

"""

from __future__ import absolute_import
import pyclaw.limiters.weno.reconstruct as recon
from six.moves import range

def weno(k, q):
    r"""Return the *k* order WENO based reconstruction of *q*.

    The reconstruction is component based.
    """

    import numpy as np
    # XXX: this should really by a class so that the workspaces (sigma
    # and weights) can be pre-allocated and the 'getattr's can be done
    # once instead of every call

    if (k % 2) == 0:
        raise ValueError('even order WENO reconstructions are not supported')

    k = (k+1)/2
    sigma = np.zeros((q.shape[1], k))
    weights = np.zeros((q.shape[1], k))

    ql = np.zeros(q.shape)
    qr = np.zeros(q.shape)

    try:
        smoothness    = getattr(recon, 'smoothness_k' + str(k))
        weights_l     = getattr(recon, 'weights_left_k' + str(k))
        weights_r     = getattr(recon, 'weights_right_k' + str(k))
        reconstruct_l = getattr(recon, 'reconstruct_left_k' + str(k))
        reconstruct_r = getattr(recon, 'reconstruct_right_k' + str(k))
    except:
        raise ValueError('%d order WENO reconstructions are not supported' % (2*k-1))


    for m in range(q.shape[0]):
        smoothness(q[m,:], sigma)

        weights_l(sigma, weights)
        reconstruct_l(q[m,:], weights, ql[m,:])

        weights_r(sigma, weights)
        reconstruct_r(q[m,:], weights, qr[m,:])


    # XXX: copy ghost-cells.  i'm not sure why this is necessary, but
    # it make the acoustics examples in the implicit time-stepping
    # branch work properly.

    ql[:,:k-1]  = ql[:,-2*k+2:-k+1]
    ql[:,-k+1:] = ql[:,k-1:2*k-2]

    qr[:,:k-1]  = qr[:,-2*k+2:-k+1]
    qr[:,-k+1:] = qr[:,k-1:2*k-2]

    return ql, qr
