#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing petclaw grid.

:Authors:
    Amal Alghamdi
    David Ketcheson
    Aron Ahmadia
"""

import pyclaw.grid
# We don't use Grid directly but we need it to appear in this namespace:
from pyclaw.grid import Grid

# ============================================================================
#  Dimension Object
# ============================================================================
class Dimension(pyclaw.grid.Dimension):
    r"""
    Basic class representing a dimension of a Grid object

    The only difference between PyClaw and PetClaw grids are the
    boundary conditions.
    
    :Initialization:
    
    Input:
     - *name* - (string) string Name of dimension
     - *lower* - (float) Lower extent of dimension
     - *upper* - (float) Upper extent of dimension
     - *n* - (int) Number of grid cells
     - *units* - (string) Type of units, used for informational purposes only
        
    Output:
     - (:class:`Dimension`) - Initialized Dimension object

    We could use ng, nstart, and nend in the edge and center properties in 
    PyClaw and then nothing would need to be overridden here.  But it would
    put some parallel-ish code in PyClaw, which is undesirable.
    """
    @property
    def ng(self):
        r"""Size of this processes' piece of grid in given dimension."""
        return self.nend-self.nstart

    @property
    def edge(self):
        r"""(ndarrary(:)) - Location of all grid cell edge coordinates
        for this dimension"""
        import numpy as np
        if self._edge is None:
            self._edge = np.empty(self.ng+1)
            for i in xrange(self.nstart,self.nend+1):
                self._edge[i] = self.lower + i*self.d
        return self._edge
    _edge = None

    @property
    def center(self):
        r"""(ndarrary(:)) - Location of all grid cell center coordinates
        for this dimension"""
        import numpy as np
        if self._center is None:
            self._center = np.empty(self.ng)
            for i in xrange(self.nstart,self.nend):
                self._center[i-self.nstart] = self.lower + (i+0.5)*self.d
        return self._center
    _center = None
