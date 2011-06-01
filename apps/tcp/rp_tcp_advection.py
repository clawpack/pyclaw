#!/usr/bin/env python
# encoding: utf-8
r"""
Variable coefficient 1D advection Riemann solver
for a system of 2 equations, used for modeling TCP Vegas

:Authors:
    David I. Ketcheson (Feb. 2011)
"""

import numpy as np

def rp_tcp_advection_1d(q_l,q_r,aux_l,aux_r,aux_global):
    r"""Basic 1d 2-equation advection riemann solver
    """
    
    # Riemann solver constants
    meqn = 2
    mwaves = 2
    
    # Number of Riemann problems we are solving
    nrp = q_l.shape[1]
    
    # Return values
    wave = np.empty( (meqn,mwaves,nrp) )
    s = np.empty( (mwaves,nrp) )
    amdq = np.zeros( (meqn,nrp) )
    apdq = np.zeros( (meqn,nrp) )
    
    wave[0,0,:] = q_r[0,:] - q_l[0,:]
    wave[1,0,:] = 0.
    wave[0,1,:] = 0.
    wave[1,1,:] = q_r[1,:] - q_l[1,:]

    s[0,:] = aux_l[0,:]
    s[1,:] = aux_l[1,:]

    apdq[0,:] = s[0,:] * wave[0,0,:]
    apdq[1,:] = s[0,:] * wave[1,0,:]
    amdq[0,:] = s[1,:] * wave[0,1,:]
    amdq[1,:] = s[1,:] * wave[1,1,:]

    return wave, s, amdq, apdq
