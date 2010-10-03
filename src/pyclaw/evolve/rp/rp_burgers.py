#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for Burgers equation

.. math::
    u_t + \left ( \frac{1}{2} u^2 \right)_x = 0

:Authors:
    Kyle T. Mandli (2009-2-4): Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

def rp_burgers_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""
    Riemann solver for Burgers equation in 1d
         
    *aux_global* should contain -
     - *efix* - (bool) Whether a entropy fix should be used, if not present, 
       false is assumed
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2009-2-4)
    """
        
    nrp = q_l.shape[0]
    meqn = 1
    mwaves = 1
    
    # Output arrays
    wave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.empty( (nrp, meqn) )
    apdq = np.empty( (nrp, meqn) )

    # Basic solve
    wave[:,:,0] = q_r - q_l
    s[:,0] = 0.5 * (q_r[:,0] + q_l[:,0])

    s_index = np.zeros((nrp,2))
    s_index[:,0] = s[:,0]
    amdq[:,0] = np.min(s_index,axis=1) * wave[:,0,0]
    apdq[:,0] = np.max(s_index,axis=1) * wave[:,0,0]
    
    # Compute entropy fix
    if aux_global['efix']:
        transonic = (q_l[:,0] < 0.0) * (q_r[:,0] > 0.0)
        amdq[transonic,0] = -0.5 * q_l[transonic,0]**2
        apdq[transonic,0] = 0.5 * q_r[transonic,0]**2

    return wave, s, amdq, apdq
