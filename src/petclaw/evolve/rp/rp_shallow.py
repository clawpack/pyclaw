#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for the shallow water equations.
    
The available solvers are:
 * Roe - Use Roe averages to caluclate the solution to the Riemann problem
 * HLL - Use a HLL solver
 * Exact - Use a newton iteration to calculate the exact solution to the 
        Riemann problem

.. math:: 
    q_t + f(q)_x = 0
  
where 

.. math:: 
    q(x,t) = \left [ \begin{array}{c} h \\ h u \end{array} \right ],
  
the flux function is 

.. math:: 
    f(q) = \left [ \begin{array}{c} h u \\ hu^2 + 1/2 g h^2 \end{array}\right ].

and :math:`h` is the water column height, :math:`u` the velocity and :math:`g`
is the gravitational acceleration.

:Authors:
    Kyle T. Mandli (2009-02-05): Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

def rp_shallow_roe_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""
    Roe shallow water solver in 1d::
    
        ubar = (sqrt(u_l) + sqrt(u_r)) / (sqrt(h_l) + sqrt(h_r))
        cbar = sqrt( 0.5 * g * (h_l + h_r))
         
        W_1 = |      1      |  s_1 = ubar - cbar
              | ubar - cbar |
         
        W_2 = |      1      |  s_1 = ubar + cbar
              | ubar + cbar |
          
        a1 = 0.5 * ( - delta_hu + (ubar + cbar) * delta_h ) / cbar
        a2 = 0.5 * (   delta_hu - (ubar - cbar) * delta_h ) / cbar
    
    *aux_global* should contain:
     - *g* - (float) Gravitational constant
     - *efix* - (bool) Boolean as to whether a entropy fix should be used, if 
       not present, false is assumed
            
    :Version: 1.0 (2009-02-05)
    """
    
    # Array shapes
    nrp = q_l.shape[0]
    meqn = 2
    mwaves = 2
    
    # Output arrays
    wave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.zeros( (nrp, meqn) )
    apdq = np.zeros( (nrp, meqn) )

    # Compute roe-averaged quantities
    ubar = ( (q_l[:,1]/np.sqrt(q_l[:,0]) + q_r[:,1]/np.sqrt(q_r[:,0])) /
        (np.sqrt(q_l[:,0]) + np.sqrt(q_r[:,0])) )
    cbar = np.sqrt(0.5 * aux_global['g'] * (q_l[:,0] + q_r[:,0]))
        
    # Compute Flux structure    
    delta = q_r - q_l
    a1 = 0.5 * (-delta[:,1] + (ubar + cbar) * delta[:,0]) / cbar
    a2 = 0.5 * ( delta[:,1] - (ubar - cbar) * delta[:,0]) / cbar
        
    # Compute each family of waves
    wave[:,0,0] = a1
    wave[:,1,0] = a1 * (ubar - cbar)
    s[:,0] = ubar - cbar
    
    wave[:,0,1] = a2
    wave[:,1,1] = a2 * (ubar + cbar)
    s[:,1] = ubar + cbar
        
    if aux_global['efix']:
        raise NotImplementedError("Entropy fix has not been implemented.")
    else:
        s_index = np.zeros((nrp,2))
        for m in xrange(meqn):
            for mw in xrange(mwaves):
                s_index[:,0] = s[:,mw]
                amdq[:,m] += np.min(s_index,axis=1) * wave[:,m,mw]
                apdq[:,m] += np.max(s_index,axis=1) * wave[:,m,mw]
            
    return wave, s, amdq, apdq
    
def rp_shallow_hll_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""
    HLL shallow water solver ::
    
         
        W_1 = Q_hat - Q_l    s_1 = min(u_l-c_l,u_l+c_l,lambda_roe_1,lambda_roe_2)
        W_2 = Q_r - Q_hat    s_2 = max(u_r-c_r,u_r+c_r,lambda_roe_1,lambda_roe_2)
    
        Q_hat = ( f(q_r) - f(q_l) - s_2 * q_r + s_1 * q_l ) / (s_1 - s_2)
    
    *aux_global* should contain:
     - *g* - (float) Gravitational constant
            
    :Version: 1.0 (2009-02-05)
    """
    
    # Array shapes
    nrp = q_l.shape[0]
    meqn = 2
    mwaves = 2
    
    # Output arrays
    wave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.zeros( (nrp, meqn) )
    apdq = np.zeros( (nrp, meqn) )

    # Compute Roe and right and left speeds
    ubar = ( (q_l[:,1]/np.sqrt(q_l[:,0]) + q_r[:,1]/np.sqrt(q_r[:,0])) /
        (np.sqrt(q_l[:,0]) + np.sqrt(q_r[:,0])) )
    cbar = np.sqrt(0.5 * aux_global['g'] * (q_l[:,0] + q_r[:,0]))
    u_r = q_r[:,1] / q_r[:,0]
    c_r = np.sqrt(aux_global['g'] * q_r[:,0])
    u_l = q_l[:,1] / q_l[:,0]
    c_l = np.sqrt(aux_global['g'] * q_l[:,0])

    # Compute Einfeldt speeds
    s_index = np.empty((nrp,4))
    s_index[:,0] = ubar+cbar
    s_index[:,1] = ubar-cbar
    s_index[:,2] = u_l + c_l
    s_index[:,3] = u_l - c_l
    s[:,0] = np.min(s_index,axis=1)
    s_index[:,2] = u_r + c_r
    s_index[:,3] = u_r - c_r
    s[:,1] = np.max(s_index,axis=1)

    # Compute middle state
    q_hat = np.empty((nrp,2))
    q_hat[:,0] = ((q_r[:,1] - q_l[:,1] - s[:,1] * q_r[:,0] 
                            + s[:,0] * q_l[:,0]) / (s[:,0] - s[:,1]))
    q_hat[:,1] = ((q_r[:,1]**2/q_r[:,0] + 0.5 * aux_global['g'] * q_r[:,0]**2
                - (q_l[:,1]**2/q_l[:,0] + 0.5 * aux_global['g'] * q_l[:,0]**2)
                - s[:,1] * q_r[:,1] + s[:,0] * q_l[:,1]) / (s[:,0] - s[:,1]))

    # Compute each family of waves
    wave[:,:,0] = q_hat - q_l
    wave[:,:,1] = q_r - q_hat
    
    # Compute variations
    s_index = np.zeros((nrp,2))
    for m in xrange(meqn):
        for mw in xrange(mwaves):
            s_index[:,0] = s[:,mw]
            amdq[:,m] += np.min(s_index,axis=1) * wave[:,m,mw]
            apdq[:,m] += np.max(s_index,axis=1) * wave[:,m,mw]
            
    return wave, s, amdq, apdq
    
def rp_shallow_exact_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""
    Exact shallow water Riemann solver
    
    .. warning::
        This solver has not been implemented.
    
    """
    raise NotImplementedError("The exact swe solver has not been implemented.")
