#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing petclaw.geometry.

:Authors:
    Amal Alghamdi
    David Ketcheson
    Aron Ahmadia
"""

import pyclaw.geometry
# We don't use Patch directly but we need it to appear in this namespace:
from pyclaw.geometry import Patch

# ============================================================================
#  Dimension Object
# ============================================================================
class Dimension(pyclaw.geometry.Dimension):
    r"""
    Basic class representing a dimension of a Patch object

    The only difference between PyClaw and PetClaw patchs are the
    boundary conditions.
    
    :Initialization:
    
    Input:
     - *name* - (string) string Name of dimension
     - *lower* - (float) Lower extent of dimension
     - *upper* - (float) Upper extent of dimension
     - *n* - (int) Number of patch cells
     - *units* - (string) Type of units, used for informational purposes only
        
    Output:
     - (:class:`Dimension`) - Initialized Dimension object

    We could use ng, nstart, and nend in the edge and center properties in 
    PyClaw and then nothing would need to be overridden here.  But it would
    put some parallel-ish code in PyClaw, which is undesirable.
    """
    @property
    def ng(self):
        r"""Size of this processes' piece of patch in given dimension."""
        return self.nend-self.nstart

    @property
    def edges(self):
        r"""(ndarrary(:)) - Location of all patch cell edge coordinates
        for this dimension"""
        import numpy as np
        if self._edges is None:
            self._edges = np.empty(self.ng+1)
            for i in xrange(self.nstart,self.nend+1):
                self._edges[i] = self.lower + i*self.delta
        return self._edges
    _edges = None

    @property
    def centers(self):
        r"""(ndarrary(:)) - Location of all patch cell center coordinates
        for this dimension"""
        import numpy as np
        if self._centers is None:
            self._centers = np.empty(self.ng)
            for i in xrange(self.nstart,self.nend):
                self._centers[i-self.nstart] = self.lower + (i+0.5)*self.delta
        return self._centers
    _centers = None
