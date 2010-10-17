#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for constant coefficient acoustics.

.. math:: q_t + A q_x = 0
  
where 

.. math:: 
    q(x,t) = \left [ \begin{array}{c} p(x,t) \\ u(x,t) \end{array} \right ]
  
and the coefficient matrix is 

.. math::

    A = \left [\begin{matrix}
    0      & K\\
    1/\rho & 0 
    \end{matrix} \right ]

The parameters :math:`\rho =` density and :math:`K =` bulk modulus are used
to calculate the impedence :math:`= Z` and speed of sound `= c`.

:Authors:
    Kyle T. Mandli (2009-02-03): Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

def rp_acoustics_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""
    Basic 1d acoustics riemann solver
    
    *aux_global* is expected to contain -
     - *zz* - (float) Impedence
     - *cc* - (float) Speed of sound
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2009-02-03)
    """
    
    meqn = 2
    mwaves = 2
    # Convenience
    nrp = np.size(q_l,0)

    # Return values
    wave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.empty( (nrp, meqn) )
    apdq = np.empty( (nrp, meqn) )
    
    # Local values
    delta = np.empty(np.shape(q_l))
    
    delta = q_r - q_l
    a1 = (-delta[:,0] + aux_global['zz']*delta[:,1]) / (2.0 * aux_global['zz'])
    a2 = (delta[:,0] + aux_global['zz']*delta[:,1]) / (2.0 * aux_global['zz'])
        
    # Compute the waves
    # 1-Wave
    wave[:,0,0] = -a1 * aux_global['zz']
    wave[:,1,0] = a1
    s[:,0] = -aux_global['cc']
        
    # 2-Wave
    wave[:,0,1] = a2 * aux_global['zz']
    wave[:,1,1] = a2
    s[:,1] = aux_global['cc']
    
    # Compute the left going and right going fluctuations
    for m in xrange(meqn):
        amdq[:,m] = s[:,0] * wave[:,m,0]
        apdq[:,m] = s[:,1] * wave[:,m,1]
    
    return wave, s, amdq, apdq
