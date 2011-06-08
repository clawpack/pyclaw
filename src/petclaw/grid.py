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

import numpy as np

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
    put some parallel code in PyClaw.
    """
    def ng():
        doc = r"""Size of this processes' piece of grid in given dimension."""
        def fget(self):
            return self.nend-self.nstart
        return locals()
    ng = property(**ng())

    def edge():
        doc = r"""(ndarrary(:)) - Location of all grid cell edge coordinates
        for this dimension"""
        def fget(self): 
            import numpy as np
            if self._edge is None:
                self._edge = np.empty(self.ng+1)
                for i in xrange(self.nstart,self.nend+1):
                    self._edge[i] = self.lower + i*self.d
            return self._edge
        return locals()
    edge = property(**edge())
    _edge = None

    def center():
        doc = r"""(ndarrary(:)) - Location of all grid cell center coordinates
        for this dimension"""
        def fget(self): 
            import numpy as np
            if self._center is None:
                self._center = np.empty(self.ng)
                for i in xrange(self.nstart,self.nend):
                    self._center[i-self.nstart] = self.lower + (i+0.5)*self.d
            return self._center
        return locals()
    center = property(**center())
    _center = None



# ============================================================================
#  petclaw Grid object definition
# ============================================================================
class Grid(pyclaw.grid.Grid):
    r"""
    Basic representation of a single grid in petclaw

    The only difference between PetClaw grid and PyClaw grid is
    the definition of q(), local_n(), __getstate__(), 
    and __setstate__().
    
    :Dimension information:
    
        Each dimension has an associated name with it that can be accessed via
        that name such as ``grid.x.n`` which would access the x dimension's
        number of grid cells.
    
    :Global Grid information:
    
        Each grid has a value for :attr:`level` and :attr:`gridno`.
        
    :Grid Data:
    
        The array :attr:`capa` has variable 
        extents based on the set of dimensions present and the values of 
        :attr:`meqn` and :attr:`maux`.  
        The :attr:`capa` array is 
        initially set to all ``1.0`` and needs to be manually set.
        
    :Properties:

        If the requested property has multiple values, a list will be returned
        with the corresponding property belonging to the dimensions in order.
         
    :Initialization:
    
        Input:
         - *dimensions* - (list of :class:`Dimension`) Dimensions that are to 
           be associated with this grid
            
        Output:
         - (:class:`Grid`) Initialized grid object
    """
