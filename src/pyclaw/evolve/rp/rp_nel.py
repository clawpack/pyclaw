#!/usr/bin/env python
# encoding: utf-8
r"""
F-wave Riemann solver for nonlinear elasticity in heterogeneous media

.. math:: q_t + f(q,x)_x = 0
  
where 

.. math:: 
    q(x,t) = \left [ \begin{array}{c} \epsilon(x,t) \\ \rho(x) u(x,t) \end{array} \right ]
  
and the flux vector is

.. math::

    f(q,x) = \left [ \begin{array}{c} -u \\ \sigma(\epsilon,x) \end{array} \right ]

:Authors:
    David I. Ketcheson (2010-11-06): Initial version
"""
# ============================================================================
#      Copyright (C) 2010 David I. Ketcheson <david.ketcheson@kaust.edu.sa>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

def rp_nel_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""
    1d nonlinear elasticity riemann solver
    
    *aux* is expected to contain -
     - aux[i,0] - density in cell i
     - aux[i,1] - bulk modulus in cell i
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2010-11-06)
    """
    
    meqn = 2
    mwaves = 2
    # Convenience
    nrp = np.size(q_l,0)

    # Set up arrays for return values
    fwave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.empty( (nrp, meqn) )
    apdq = np.empty( (nrp, meqn) )
    
    #Linearized bulk modulus, sound speed, and impedance:
    bulkl = sigmap(q_l[:,0],aux_l[:,1])
    bulkr = sigmap(q_r[:,0],aux_r[:,1])
    cl = np.sqrt(bulkl/aux_l[:,0])
    cr = np.sqrt(bulkr/aux_r[:,0])
    zl = cl*aux_l[:,0]
    zr = cr*aux_r[:,0]

    #Jumps:
    du   = q_r[:,1]/aux_r[:,0]-q_l[:,1]/aux_l[:,0]
    dsig = sigma(q_r[:,0],aux_r[:,1]) - sigma(q_l[:,0],aux_l[:,1])

    b1 = - (zr*du + dsig) / (zr+zl)
    b2 = - (zl*du - dsig) / (zr+zl)

    # Compute the f-waves
    # 1-Wave
    fwave[:,0,0] = b1
    fwave[:,1,0] = b1 * zl
    s[:,0] = -cl
        
    # 2-Wave
    fwave[:,0,1] = b2
    fwave[:,1,1] = b2 * (-zr)
    s[:,1] = cr
    
    # Compute the left going and right going fluctuations
    for m in xrange(meqn):
        amdq[:,m] = fwave[:,m,0]
        apdq[:,m] = fwave[:,m,1]
    
    return fwave, s, amdq, apdq

def sigma(eps,K):
    return np.exp(K*eps)-1.0

def sigmap(eps,K):
    return K*np.exp(K*eps)
