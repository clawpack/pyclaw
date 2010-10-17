#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for the Euler equations

This module contains Riemann solvers for the Euler equations which have the
form (in 1d):

.. math:: 
    q_t + f(q)_x = 0
  
where 

.. math:: 
    q(x,t) = \left [ \begin{array}{c} \rho \\ \rho u \\ E \end{array} \right ],
  
the flux function is 

.. math:: 
    f(q) = \left [ \begin{array}{c} \rho u \\ \rho u^2 + p \\ u(E+p) \end{array}\right ].

and :math:`\rho` is the density, :math:`u` the velocity, :math:`E` is the 
energy and :math:`p` is the pressure.

Unless otherwise noted, the ideal gas equation of state is used:

.. math::
    E = (\gamma - 1) \left (E - \frac{1}{2}\rho u^2 \right)

:Authors:
    Kyle T. Mandli (2009-06-26): Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

def rp_euler_roe_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""
    Roe Euler solver in 1d
    
    *aug_global* should contain -
     - *gamma* - (float) Ratio of the heat capacities
     - *gamma1* - (float) :math:`1 - \gamma`
     - *efix* - (bool) Whether to use an entropy fix or not
    
    See :ref:`pyclaw_rp` for more details.
    
    :Version: 1.0 (2009-6-26)
    """
    
    # Problem dimensions
    meqn = 3
    mwaves = 3
    nrp = np.size(q_l,0)

    # Return values
    wave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.zeros( (nrp, meqn) )
    apdq = np.zeros( (nrp, meqn) )
    
    # Solver parameters
    gamma = aux_global['gamma']
    gamma1 = aux_global['gamma1']

    # Calculate Roe averages
    rhsqrtl = np.sqrt(q_l[...,0])
    rhsqrtr = np.sqrt(q_r[...,0])
    pl = gamma1 * (q_l[...,2] - 0.5 * (q_l[...,1]**2) / q_l[...,0])
    pr = gamma1 * (q_r[...,2] - 0.5 * (q_r[...,1]**2) / q_r[...,0])
    rhsq2 = rhsqrtl + rhsqrtr
    u = q_l[...,1] / rhsqrtl + q_r[...,1] / rhsqrtr
    enthalpy = ((q_l[...,2] + pl) / rhsqrtl + (q_r[...,2] + pr) / rhsqrtr) / rhsq2
    a = np.sqrt(gamma1 * (enthalpy - 0.5 * u**2))

    # Find eigenvector coefficients
    delta = q_r - q_l
    a2 = gamma1 / a**2 * ((enthalpy -u**2)*delta[...,0] + u*delta[...,1] - delta[...,2])
    a3 = (delta[...,1] + (a-u) * delta[...,0] - a*a2) / (2.0*a)
    a1 = delta[...,0] - a2 - a3
    
    # Compute the waves
    wave[...,0,0] = a1
    wave[...,1,0] = a1 * (u-a)
    wave[...,2,0] = a1 * (enthalpy - u*a)
    s[...,0] = u - a
    
    wave[...,0,1] = a2
    wave[...,1,1] = a2 * u
    wave[...,2,1] = a2 * 0.5 * u**2
    s[...,1] = u
    
    wave[...,0,2] = a3
    wave[...,1,2] = a3 * (u+a)
    wave[...,2,2] = a3 * (enthalpy + u*a)
    s[...,2] = u + a
    
    # Entropy fix
    if aux_global['efix']:
        raise NotImplementedError("Entropy fix has not been implemented!")
    else:
        # Godunov update
        s_index = np.zeros((nrp,2))
        for m in xrange(meqn):
            for mw in xrange(mwaves):
                s_index[:,0] = s[:,mw]
                amdq[:,m] += np.min(s_index,axis=1) * wave[:,m,mw]
                apdq[:,m] += np.max(s_index,axis=1) * wave[:,m,mw]
    

    return wave,s,amdq,apdq