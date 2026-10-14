#!/usr/bin/env python
# encoding: utf-8
r"""
Library of limiter functions to be applied to waves

This module contains all of the standard limiters found in clawpack.  To
use any of the limiters, use the function limit to limit the appropriate
waves.  Refer to each limiter and the function limit's doc strings.

This is a list of the provided limiters and their corresponding method number,
note that some of the limiters actually correspond to a more general function
which can be controlled more directly.  Refer to the limiter function and its
corresponding documentation for details.

CFL Independent Limiters
''''''''''''''''''''''''
1. minmod - :func:`minmod_limiter`
2. superbee - :func:`superbee_limiter`
3. van leer - :math:`(r + |r|) / (1 + |r|)`
4. mc - :func:`mc_limiter`
5. Beam-warming - :math:`r`
6. Frommm - :math:`1/2 (1 + r)`
7. Albada 2 - :math:`(r^2 + r) / (1 + r^2)`
8. Albada 3 - :math:`1/2 (1+r) (1 - (|1-r|^3) / (1+|r|^3))`
9. van Leer with Klein sharpening, k=2 - 
   :func:`van_leer_klein_sharpening_limiter`

CFL Dependent Limiters
''''''''''''''''''''''
10. Roe's linear third order scheme - :math:`1 + (r-1) (1 + cfl) / 3`
11. Arora-Roe (= limited version of the linear third order scheme) - 
    :func:`arora_roe`
12. Theta Limiter, theta=0.95 (safety on nonlinear waves) - 
    :func:`theta_limiter`
13. Theta Limiter, theta=0.75 - :func:`theta_limiter`
14. Theta Limiter, theta=0.5 - :func:`theta_limiter`
15. CFL-Superbee (Roe's Ultrabee)  - :func:`cfl_superbee`
16. CFL-Superbee (Roe's Ultrabee) with theta=0.95 (nonlinear waves) - 
    :func:`cfl_superbee_theta`
17. beta=2/3 limiter - :func:`beta_limiter`
18. beta=2/3 limiter with theta=0.95 (nonlinear waves) - :func:`beta_limiter`
19. Hyperbee - :func:`hyperbee_limiter`
20. SuperPower - :func:`superpower_limiter`
21. Cada-Torrilhon modified - :func:`cada_torrilhon_limiter`
22. Cada-Torrilhon modified, version for nonlinear waves - 
    :func:`cada_torrilhon_limiter_nonlinear`
23. upper bound limiter (1st order) - :func:`upper_bound_limiter`
    
All limiters have the same function call signature:
    :Input:
     - *r* - (ndarray(:)) 
     - *cfl* - (ndarray(:)) Local CFL number
     
    :Output:
     - (ndarray(:)) - 

Newer limiters are based on work done by Friedemann Kemm [kemm_2009]_, paper 
in review.
    
:Authors:
    Kyle Mandli and Randy LeVeque (2008-08-21) Initial version
    
    Kyle Mandli (2009-07-05) Added CFL depdendent limiters
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#      Copyright (C) 2009 Randall J. LeVeque <rjl@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

minmod = 1
superbee = 2
vanleer = 3
MC = 4

import numpy as np

def limit(num_eqn,wave,s,limiter,dtdx):
    r"""
    Apply a limiter to the waves

    Function that limits the given waves using the methods contained
    in limiter.  This is the vectorized version of the function acting on a 
    row of waves at a time.
    
    :Input:
     - *wave* - (ndarray(:,num_eqn,num_waves)) The waves at each interface
     - *s* - (ndarray(:,num_waves)) Speeds for each wave
     - *limiter* - (``int`` list) Array of type ``int`` determining which 
         limiter to use
     - *dtdx* - (ndarray(:)) :math:`\Delta t / \Delta x` ratio, used for CFL 
        dependent limiters
        
    :Output:
     - (ndarray(:,num_eqn,num_waves)) - Returns the limited waves

    :Version: 1.1 (2009-07-05)
    """
    
    # wave_norm2 is the sum of the squares along the num_eqn axis,
    # so the norm of the cell i for wave number j is addressed 
    # as wave_norm2[i,j]
    wave_norm2 = np.sum(np.square(wave),axis=0)
    wave_zero_mask = np.array((wave_norm2 == 0), dtype=float)
    wave_nonzero_mask = (1.0-wave_zero_mask)

    # dotls contains the products of adjacent cell values summed
    # along the num_eqn axis.  For reference, dotls[0,:,:] is the dot
    # product of the 0 cell and the 1 cell.
    dotls = np.sum(wave[:,:,1:]*wave[:,:,:-1],axis=0)

    # array containing ones where s > 0, zeros elsewhere
    spos = np.array(s > 0.0, dtype=float)[:,1:-1]

    # Here we construct a masked array, then fill the empty values with 0,
    # this is done in case wave_norm2 is 0 or close to it
    # Take upwind dot product
    r = np.ma.array((spos*dotls[:,:-1] + (1-spos)*dotls[:,1:]))
    # Divide it by the norm**2
    r /= np.ma.array(wave_norm2[:,1:-1])
    # Fill the rest of the array
    r.fill_value = 0
    r = r.filled()
    
    for mw in range(wave.shape[1]):
        # skip waves that are marked as not needing a limiter
        limit_func = limiter_functions.get(limiter[mw])
        if limit_func is not None:
            for m in range(num_eqn):
                cfl = np.abs(s[mw,1:-1]*(dtdx[1:-2]*spos[mw,:] 
                                        + (1-spos[mw,:])*dtdx[2:-1]))
                wlimitr = limit_func(r[mw,:],cfl)
                wave[m,mw,1:-1] = wave[m,mw,1:-1]*wave_zero_mask[mw,1:-1] \
                    + wlimitr * wave[m,mw,1:-1] * wave_nonzero_mask[mw,1:-1]

    return wave

def minmod_limiter(r,cfl):
    r"""
    Minmod vectorized limiter
    """
    a = np.ones((2,len(r)))
    b = np.zeros((2,len(r)))
    
    a[1,:] = r
    
    b[1,:] = np.minimum(a[0],a[1])
    
    return np.maximum(b[0],b[1])
    

def superbee_limiter(r,cfl):
    r"""
    Superbee vectorized limiter
    """
    a = np.ones((2,len(r)))
    b = np.zeros((2,len(r)))
    c = np.zeros((3,len(r)))

    a[1,:] = 2.0*r
    
    b[1,:] = r

    c[1,:] = np.minimum(a[0], a[1])
    c[2,:] = np.minimum(b[0], b[1])

    return np.max(c,axis=0)
    
def mc_limiter(r,cfl):
    r"""
    MC vectorized limiter
    """
    a = np.empty((3,len(r)))
    b = np.zeros((2,len(r)))

    a[0,:] = (1.0 + r) / 2.0
    a[1,:] = 2
    a[2,:] = 2.0 * r
    
    b[1,:] = np.min(a,axis=0)

    return np.maximum(b[0], b[1])
    
def van_leer_klein_sharpening_limiter(r,cfl):
    r"""
    van Leer with Klein sharpening, k=2
    """
    a = np.ones((2,len(r))) * 1.e-5
    a[0,:] = r
    
    rcorr = np.maximum(a[0],a[1])
    a[1,:] = 1/rcorr
    sharg = np.minimum(a[0],a[1])
    sharp = 1.0 + sharg * (1.0 - sharg) * (1.0 - sharg**2)
    
    return (r+np.abs(r)) / (1.0 + np.abs(r)) * sharp

def arora_roe(r,cfl):
    r"""
    Arora-Roe limiter, limited version of the linear third order scheme
    """
    caut = 0.99
    
    a = np.empty((3,len(r)))
    b = np.zeros((2,len(r)))
    
    s1 = (caut * 2.0 / cfl)
    s2 = (1.0 + cfl) / 3.0
    phimax = caut * 2.0 / (1.0 - cfl)
    
    a[0,:] = s1 * r
    a[1,:] = 1.0 + s2 * (r - 1.0)
    a[2,:] = phimax
    b[1,:] = np.min(a,axis=0)
    
    return np.maximum(b[0],b[1])

def theta_limiter(r,cfl,theta=0.95):
    r"""
    Theta limiter
    
    Additional Input:
     - *theta* =
    """
    a = np.empty((2,len(r)))
    b = np.empty((3,len(r)))
    
    a[0,:] = 0.001
    a[1,:] = cfl
    cfmod1 = np.max(a,axis=0)
    a[0,:] = 0.999
    cfmod2 = np.min(a,axis=0)
    s1 = 2.0 / cfmod1
    s2 = (1.0 + cfl) / 3.0
    phimax = 2.0 / (1.0 - cfmod2)
    
    a[0,:] = (1.0 - theta) * s1
    a[1,:] = 1.0 + s2 * (r - 1.0)
    left = np.maximum(a[0],a[1])
    a[0,:] = (1.0 - theta) * phimax * r
    a[1,:] = theta * s1 * r
    middle = np.maximum(a[0],a[1])
    
    b[0,:] = left
    b[1,:] = middle
    b[2,:] = theta*phimax
    
    return np.min(b,axis=0)
    
def cfl_superbee(r,cfl):
    r"""
    CFL-Superbee (Roe's Ultrabee) without theta parameter
    """
    a = np.empty((2,len(r)))
    b = np.zeros((3,len(r)))
    
    a[0,:] = 0.001
    a[1,:] = cfl
    cfmod1 = np.maximum(a[0],a[1])
    a[0,:] = 0.999
    cfmod2 = np.minimum(a[0],a[1])
    
    a[0,:] = 1.0
    a[1,:] = 2.0 * r / cfmod1
    b[1,:] = np.minimum(a[0],a[1])
    a[0,:] = 2.0/(1-cfmod2)
    a[1,:] = r
    b[2,:] = np.minimum(a[0],a[1])
    
    return np.max(b,axis=0)

    
def cfl_superbee_theta(r,cfl,theta=0.95):
    r"""
    CFL-Superbee (Roe's Ultrabee) with theta parameter
    """
    a = np.empty((2,len(r)))
    b = np.zeros((2,len(r)))
    
    a[0,:] = 0.001
    a[1,:] = cfl
    cfmod1 = np.maximum(a[0],a[1])
    a[0,:] = 0.999
    cfmod2 = np.minimum(a[0],a[1])

    s1 = theta * 2.0 / cfmod1
    phimax = theta * 2.0 / (1.0 - cfmod2)

    a[0,:] = s1*r
    a[1,:] = phimax
    b[1,:] = np.minimum(a[0],a[1])
    ultra = np.maximum(b[0],b[1])
    
    a[0,:] = ultra
    b[0,:] = 1.0
    b[1,:] = r
    a[1,:] = np.maximum(b[0],b[1])
    return np.min(a,axis=0)

def beta_limiter(r,cfl,theta=0.95,beta=0.66666666666666666):
    r"""
    Modification of CFL Superbee limiter with theta and beta parameters
    
    Additional Input:
     - *theta*
     - *beta*
    """
    a = np.empty((2,len(r)))
    b = np.zeros((2,len(r)))
    
    a[0,:] = 0.001
    a[1,:] = cfl
    cfmod1 = np.maximum(a[0],a[1])
    a[0,:] = 0.999
    cfmod2 = np.minimum(a[0],a[1])
    
    s1 = theta * 2.0 / cfmod1
    s2 = (1.0 + cfl) / 3.0
    phimax = theta * 2.0 / (1.0 - cfmod2)
    
    a[0,:] = s1*r
    a[1,:] = phimax
    b[1,:] = np.min(a)
    ultra = np.max(b)
    
    a[0,:] = 1.0 + (s2 - beta/2.0) * (r-1.0)
    a[1,:] = 1.0 + (s2 + beta/2.0) * (r-1.0)
    b[0,:] = ultra
    b[1,:] = np.max(a)
    a[0,:] = 0.0
    a[1,:] = np.min(b)
    
    return np.max(a)


def hyperbee_limiter(r,cfl):
    r"""Hyperbee"""
    a = np.empty((2,len(r)))
    
    a[0,:] = 0.001
    a[1,:] = cfl
    cfmod1 = np.maximum(a[0],a[1])
    a[0,:] = 0.999
    cfmod2 = np.minimum(a[0],a[1])
    
    index1 = r < 0.0
    index2 = np.abs(r-1.0) < 1.0e-6
    index3 = np.abs(index1 + index2 - 1)
    master_index = index1 * 0 + index2 * 1 + index3 * 2

    rmin = r-1.0
    rdur = r/rmin
    return np.choose(master_index,[0.0,1.0,
                                    2.0 * rdur * (cfl * rmin + 1.0 - r**cfl) 
                                    / (cfmod1 * (1.0 - cfmod2) * rmin)])
    
def superpower_limiter(r,cfl,caut=1.0):
    r"""
    SuperPower limiter
    
    Additional input:
     - *caut* = Limiter parameter
    """
    s2 = (1.0 + cfl) / 3.0
    s3 = 1.0 - s2
    
    pp = (((r<=1.0) * np.abs(2.0/cfl) * caut) * 2.0 * s3 
            + (r > 1.0) * (np.abs(2.0)/(1.0 - cfl) * caut) * 2.0 * s2)
            
    rabs = np.abs(r)
    rfrac = np.abs((1.0 - rabs) / (1.0 + rabs))
    signum = np.floor(0.5 * (1.0 + np.sign(r)))
    
    return signum * (s3 + s2 * r) * (1.0 - rfrac**pp)
    
def cada_torrilhon_limiter(r,cfl,epsilon=1.0e-3):
    r"""
    Cada-Torrilhon modified
    
    Additional Input:
     - *epsilon* = 
    """
    a = np.ones((2,len(r))) * 0.95
    b = np.empty((3,len(r)))

    a[0,:] = cfl
    cfl = np.min(a)
    a[1,:] = 0.05
    cfl = np.max(a)
    
    # Multiply all parts except b[0,:] by (1.0 - epsilon) as well
    b[0,:] = 1.0 + (1+cfl) / 3.0 * (r - 1)
    b[1,:] = 2.0 * np.abs(r) / (cfl + epsilon)
    b[2,:] = (8.0 - 2.0 * cfl) / (np.abs(r) * (cfl - 1.0 - epsilon)**2)
    b[1,::2] *= (1.0 - epsilon)
    a[0,:] = np.min(b)
    a[1,:] = (-2.0 * (cfl**2 - 3.0 * cfl + 8.0) * (1.0-epsilon)
                    / (np.abs(r) * (cfl**3 - cfl**2 - cfl + 1.0 + epsilon)))
    
    return np.max(a)
    
def cada_torrilhon_limiter_nonlinear(r,cfl):
    r"""
    Cada-Torrilhon modified, version for nonlinear waves
    """
    a = np.empty((3,len(r)))
    b = np.empty((2,len(r)))
    
    s2 = (1.0 + cfl) / 3.0
    a[0,:] = 1.0 + s2 * (r - 1.0)
    a[1,:] = 2.0 * np.abs(r) / (0.6 + np.abs(r))
    a[2,:] = 5.0 / np.abs(r)
    b[0,:] = np.min(a)
    b[1,:] = -3.0 / np.abs(r)
    
    return np.max(b)
    
def upper_bound_limiter(r,cfl,theta=1.0):
    r"""
    Upper bound limiter (1st order)
    
    Additional Input:
     - *theta* =
     """
    a = np.empty((2,len(r)))
    b = np.zeros((2,len(r)))
    
    a[0,:] = 0.001
    a[1,:] = cfl
    cfmod1 = np.maximum(a[0],a[1])
    a[0,:] = 0.999
    cfmod2 = np.minimum(a[0],a[1])
    
    s1 = theta * 2.0 / cfmod1
    phimax = theta * 2.0 / (1.0 - cfmod2)
    
    a[0,:] = s1*r
    a[1,:] = phimax
    b[1,:] = np.min(a)
    
    return np.max(b)


# ============================================================================
#  Limiter function dictionary
# ============================================================================
limiter_functions = {1:minmod_limiter,
                     2:superbee_limiter,
                     3:lambda r,cfl:(r + np.abs(r)) / (1.0 + np.abs(r)),
                     4:mc_limiter,
                     5:lambda r,cfl:r,
                     6:lambda r,cfl:0.5*(1.0 + r),
                     7:lambda r,cfl:(r**2 + r) / (1.0 + r**2),
                     8:lambda r,cfl:0.5*(1+r)*(1 - np.abs(1-(r))**3/(1+np.abs(r)**3)),
                     9:van_leer_klein_sharpening_limiter,
                     10:lambda r,cfl:1.0 + (1.0 + cfl) / 3.0*(r-1.0),
                     11:arora_roe,
                     12:theta_limiter,
                     13:lambda r,cfl:theta_limiter(r,cfl,0.75),
                     14:lambda r,cfl:theta_limiter(r,cfl,0.5),
                     15:cfl_superbee,
                     16:cfl_superbee_theta,
                     17:lambda r,cfl:beta_limiter(r,cfl,theta=0.0),
                     18:beta_limiter,
                     19:hyperbee_limiter,
                     20:superpower_limiter,
                     21:cada_torrilhon_limiter,
                     22:cada_torrilhon_limiter_nonlinear,
                     23:upper_bound_limiter}
