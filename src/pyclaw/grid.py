#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing all Pyclaw solution objects

:Authors:
    Kyle T. Mandli (2008-08-07) Initial version
    David I. Ketcheson
    Amal Alghamdi
"""

import copy
import logging

import numpy as np

from data import Data

# ============================================================================
#  Default function definitions
# ============================================================================

# Default mapc2p function
def default_mapc2p(grid,x):
    r"""
    Returns the physical coordinate of the point x
    
    This is the stub function which simply returns the identity
    """
    return x


# ============================================================================
#  Dimension Object
# ============================================================================
class Dimension(object):
    r"""
    Basic class representing a dimension of a Grid object
    
    :Initialization:
    
    Input:
     - *name* - (string) string Name of dimension
     - *lower* - (float) Lower extent of dimension
     - *upper* - (float) Upper extent of dimension
     - *n* - (int) Number of grid cells
     - *units* - (string) Type of units, used for informational purposes only
       
    Output:
     - (:class:`Dimension`) - Initialized Dimension object
    """
    
    # ========== Property Definitions ========================================
    def ng():
        doc = r"""(int) - Number of grid cells in this dimension, on this process"""
        def fget(self):
            return self.n
        return locals()
    ng = property(**ng())
    def d():
        doc = r"""(float) - Size of an individual, computational grid cell"""
        def fget(self):
            return (self.upper-self.lower) / float(self.n)
        return locals()
    d = property(**d())
    def edge():
        doc = r"""(ndarrary(:)) - Location of all grid cell edge coordinates
        for this dimension"""
        def fget(self): 
            if self._edge is None:
                self._edge = np.empty(self.n+1)   
                for i in xrange(0,self.n+1):
                    self._edge[i] = self.lower + i*self.d 
            return self._edge
        return locals()
    edge = property(**edge())
    _edge = None
    def center():
        doc = r"""(ndarrary(:)) - Location of all grid cell center coordinates
        for this dimension"""
        def fget(self): 
            if self._center is None:
                self._center = np.empty(self.n)
                for i in xrange(0,self.n):
                    self._center[i] = self.lower + (i+0.5)*self.d
            return self._center
        return locals()
    center = property(**center())
    _center = None
    
    def __init__(self, *args, **kargs):
        r"""
        Creates a Dimension object
        
        See :class:`Dimension` for full documentation
        """
        
        # ========== Class Data Attributes ===================================
        self.name = 'x'
        r"""(string) Name of this coordinate dimension (e.g. 'x')"""
        self.n = None
        r"""(int) - Number of grid cells in this dimension :attr:`units`"""
        self.lower = 0.0
        r"""(float) - Lower computational grid extent"""
        self.upper = 1.0
        r"""(float) - Upper computational grid extent"""
        self.units = None
        r"""(string) Corresponding physical units of this dimension (e.g. 
        'm/s'), ``default = None``"""
        
        # Parse args
        if isinstance(args[0],float):
            self.lower = float(args[0])
            self.upper = float(args[1])
            self.n = int(args[2])
    	elif isinstance(args[0],basestring):
            self.name = args[0]
            self.lower = float(args[1])
            self.upper = float(args[2])
            self.n = int(args[3])
    	else:
    	    raise Exception("Invalid initializer for Dimension.")
        
        for (k,v) in kargs.iteritems():
            setattr(self,k,v)

        #These aren't need for PyClaw, but we set them so that
        # the PyClaw grid has the same attributes as the PetClaw
        # grid, which allows for simpler programming elsewhere.
        self.nstart = 0
        self.nend = self.n
            
        # Function attribute assignments
    

    def __str__(self):
        output = "Dimension %s" % self.name
        if self.units:
            output += " (%s)" % self.units
        output += ":  (n,d,[lower,upper]) = (%s,%s,[%s,%s])" \
            % (self.n,self.d,self.lower,self.upper)
        return output
        

# ============================================================================
#  Pyclaw Grid object definition
# ============================================================================
class Grid(object):
    r"""
    Basic representation of a single grid in Pyclaw
    
    :Dimension information:
    
        Each dimension has an associated name with it that can be accessed via
        that name such as ``grid.x.n`` which would access the x dimension's
        number of grid cells.
    
    :Global Grid information:
    
        Each grid has a value for :attr:`level` and :attr:`gridno`.
        
       
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
    
    # ========== Property Definitions ========================================
    def ndim():
        doc = r"""(int) - Number of dimensions"""
        def fget(self): return len(self._dimensions)
        return locals()
    def dimensions():
        doc = r"""(list) - List of :class:`Dimension` objects defining the 
                grid's extent and resolution"""
        def fget(self): return [getattr(self,name) for name in self._dimensions]
        return locals()
    def n():
        doc = r"""(list) - List of the total number of grid cells in each dimension"""
        def fget(self): return self.get_dim_attribute('n')
        return locals()
    def ng():
        doc = r"""(list) - List of the local (to this process)number of grid cells in each dimension"""
        def fget(self): return self.get_dim_attribute('ng')
        return locals()
    def nstart():
        doc = r"""(list) - List of the number of grid cells in each dimension"""
        def fget(self): return self.get_dim_attribute('nstart')
        return locals()
    def nend():
        doc = r"""(list) - List of the number of grid cells in each dimension"""
        def fget(self): return self.get_dim_attribute('nend')
        return locals()
    def name():
        doc = r"""(list) - List of names of each dimension"""
        def fget(self): return self._dimensions
        return locals()
    def lower():
        doc = r"""(list) - Lower coordinate extents of each dimension"""
        def fget(self): return self.get_dim_attribute('lower')
        return locals()
    def upper():
        doc = r"""(list) - Upper coordinate extends of each dimension"""
        def fget(self): return self.get_dim_attribute('upper')
        return locals()
    def d():
        doc = r"""(list) - List of computational grid cell widths"""
        def fget(self): return self.get_dim_attribute('d')
        return locals()
    def units():
        doc = r"""(list) - List of dimension units"""
        def fget(self): return self.get_dim_attribute('units')
        return locals()
    def center():
        doc = r"""(list) - List of center coordinate arrays"""
        def fget(self): return self.get_dim_attribute('center')
        return locals()
    def edge():
        doc = "List of edge coordinate arrays"
        def fget(self): return self.get_dim_attribute('edge')
        return locals()
    def p_center():
        doc = r"""(list of ndarray(...)) - List containing the arrays locating
                  the physical locations of cell centers, see 
                  :meth:`compute_p_center` for more info."""
        def fget(self):
            self.compute_p_center(self)
            return self._p_center
        return locals()
    _p_center = None
    def p_edge():
        doc = r"""(list of ndarray(...)) - List containing the arrays locating
                  the physical locations of cell edges, see 
                  :meth:`compute_p_edge` for more info."""
        def fget(self):
            self.compute_p_edge(self)
            return self._p_edge
        return locals()
    _p_edge = None
    def c_center():
        doc = r"""(list of ndarray(...)) - List containing the arrays locating
                  the computational locations of cell centers, see 
                  :meth:`compute_c_center` for more info."""
        def fget(self):
            self.compute_c_center(self)
            return self._c_center
        return locals()
    _c_center = None
    def c_edge():
        doc = r"""(list of ndarray(...)) - List containing the arrays locating
                  the computational locations of cell edges, see 
                  :meth:`compute_c_edge` for more info."""
        def fget(self):
            self.compute_c_edge(self)
            return self._c_edge
        return locals()
    _c_edge = None

       
    ndim = property(**ndim())
    dimensions = property(**dimensions())
    n = property(**n())
    ng = property(**ng())
    nstart = property(**nstart())
    nend = property(**nend())
    name = property(**name())
    lower = property(**lower())
    upper = property(**upper())
    d = property(**d())
    units = property(**units())
    center = property(**center())
    edge = property(**edge())
    p_center = property(**p_center())
    p_edge = property(**p_edge())
    c_center = property(**c_center())
    c_edge = property(**c_edge())
    
    
    # ========== Class Methods ===============================================
    def __init__(self,dimensions):
        r"""
        Instantiate a Grid object
        
        See :class:`Grid` for more info.
        """
        
        # ========== Attribute Definitions ===================================
        self.level = 1
        r"""(int) - AMR level this grid belongs to, ``default = 1``"""
        self.gridno = 1
        r"""(int) - Grid number of current grid, ``default = 0``"""
        self.mapc2p = default_mapc2p
        r"""(func) - Grid mapping function"""
        self.gauges = []
        r"""(list) - List of gauges"""
        self.gauge_files     = []
        r"""(list) - List of files to write gauge values to"""
        self.gauge_path = './_output/_gauges/'
        r"""(string) - Full path to output directory for gauges"""
        # Dimension parsing
        if isinstance(dimensions,Dimension):
            dimensions = [dimensions]
        self._dimensions = []
        for dim in dimensions:
            self.add_dimension(dim)
    
    
    def __str__(self):
        output = "Grid %s:\n" % self.gridno
        output += '\n  '.join((str(getattr(self,dim)) for dim in self._dimensions))
        output += '\n'
        return output
    
    
    def is_valid(self):
        r"""
        Checks to see if this grid is valid
        
        Nothing to do here, since both q and mbc have been moved out of grid.
            
        :Output:
         - (bool) - True if valid, false otherwise.
        
        """
        valid = True
        return valid
    
    
    # ========== Dimension Manipulation ======================================
    def add_dimension(self,dimension):
        r"""
        Add the specified dimension to this grid
        
        :Input:
         - *dimension* - (:class:`Dimension`) Dimension to be added
        """

        # Add dimension to name list and as an attribute
        self._dimensions.append(dimension.name)
        setattr(self,dimension.name,dimension)
        
        
    def get_dim_attribute(self,attr):
        r"""
        Returns a tuple of all dimensions' attribute attr
        """
        return [getattr(getattr(self,name),attr) for name in self._dimensions]
    
    
    # ========== Copy functionality ==========================================
    def __copy__(self):
        return self.__class__(self)
        
        
    def __deepcopy__(self,memo={}):
        result = self.__class__(copy.deepcopy(self.dimensions))
        result.__init__(copy.deepcopy(self.dimensions))
        
        for attr in ('level','gridno','_p_center','_p_edge',
                        '_c_center','_c_edge'):
            setattr(result,attr,copy.deepcopy(getattr(self,attr)))
        
        result.mapc2p = self.mapc2p
        
        return result
        
    
    # ========== Grid Operations =============================================
    def compute_p_center(self, recompute=False):
        r"""Calculates the :attr:`p_center` array
        
        This array is computed only when requested and then stored for later 
        use unless the recompute flag is set to True (you may want to do this
        for time dependent mappings).
        
        Access the resulting physical coordinate array via the corresponding
        dimensions or via the computational grid properties :attr:`p_center`.
        
        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        
        if recompute or not len(self._p_center) == len(self._dimensions):
            # Initialize array
            self._p_center = [None]*self.ndim

            # Special case
            if self.ndim == 1:
                self._p_center[0] = self.mapc2p(self,self.dimensions[0].center)
            # Higer dimensional calculate center arrays
            else:
                mgrid = np.lib.index_tricks.nd_grid()
                index = np.indices(self.ng)
                array_list = []
                for i,center_array in enumerate(self.get_dim_attribute('center')):
                    #We could just use indices directly and deal with
                    #numpy arrays instead of lists of numpy arrays
                    array_list.append(center_array[index[i,...]])
            
                self._p_center = self.mapc2p(self,array_list)
 

    def compute_p_edge(self, recompute=False):
        r"""Calculates the :attr:`p_edge` array
        
        This array is computed only when requested and then stored for later 
        use unless the recompute flag is set to True (you may want to do this
        for time dependent mappings).
        
        Access the resulting physical coordinate array via the corresponding
        dimensions or via the computational grid properties :attr:`p_edge`.
        
        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        
        if recompute or not len(self._p_edge) == len(self._dimensions):
            # Initialize array
            self._p_edge = [None for i in xrange(self.ndim)]

            if self.ndim == 1:        
                self._p_edge[0] = self.mapc2p(self,self.dimensions[0].edge)
            else:
                mgrid = np.lib.index_tricks.nd_grid()
                index = np.indices([n+1 for n in self.ng])
                array_list = []
                for i,edge_array in enumerate(self.get_dim_attribute('edge')):
                    #We could just use indices directly and deal with
                    #numpy arrays instead of lists of numpy arrays
                    array_list.append(edge_array[index[i,...]])
            
                self._p_edge = self.mapc2p(self,array_list)
            

    def compute_c_center(self, recompute=False):
        r"""
        Calculate the :attr:`c_center` array
        
        This array is computed only when requested and then stored for later
        use unless the recompute flag is set to True.
        
        Access the resulting computational coodinate array via the
        corresponding dimensions or via the computational grid properties
        :attr:`c_center`.
        
        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        
        if recompute or not len(self._c_center) == len(self._dimensions):
            self._c_center = [None]*self.ndim
            
            # For one dimension, the center and edge arrays are equivalent
            if self.ndim == 1:
                self._c_center[0] = self.dimensions[0].center
            else:
                # Produce ndim mesh grid function
                mgrid = np.lib.index_tricks.nd_grid()
                index = np.indices(self.ng)
                self._c_center = []
                for i,center_array in enumerate(self.get_dim_attribute('center')):
                    #We could just use indices directly and deal with
                    #numpy arrays instead of lists of numpy arrays
                    self._c_center.append(center_array[index[i,...]])

    def compute_c_edge(self, recompute=False):
        r"""
        Calculate the :attr:`c_edge` array
        
        This array is computed only when requested and then stored for later
        use unless the recompute flag is set to True.
        
        Access the resulting computational coodinate array via the
        corresponding dimensions or via the computational grid properties
        :attr:`c_edge`.
        
        :Input:
         - *recompute* - (bool) Whether to force a recompute of the arrays
        """
        if recompute or not len(self._c_edge) == len(self._dimensions):
            self._c_edge = [None for i in xrange(self.ndim)]
            
            if self.ndim == 1:
                self._c_edge[0] = self.dimensions[0].edge
            else:
                mgrid = np.lib.index_tricks.nd_grid()
                index = np.indices([n+1 for n in self.ng])
                self._c_edge = []
                for i,edge_array in enumerate(self.get_dim_attribute('edge')):
                    #We could just use indices directly and deal with
                    #numpy arrays instead of lists of numpy arrays
                    self._c_edge.append(edge_array[index[i,...]])
            

    # ========================================================================
    #  Gauges
    # ========================================================================
    def add_gauges(self,gauge_coords):
        r"""
        Determine the grid indices of each gauge and make a list of all gauges
        with their grid indices.  
        
        For PetClaw, first check whether each gauge is in the part of the grid
        corresponding to this processor.
        """
        import os
        from numpy import floor
        if not os.path.exists(self.gauge_path):
            try:
                os.makedirs(self.gauge_path)
            except OSError:
                print "gauge directory already exists, ignoring"
        
        for gauge in gauge_coords: 
            # Determine gauge locations in units of grid spacing
            gauge_ind = [int(floor(gauge[n]/self.d[n])) for n in xrange(self.ndim)]
            if all(self.nstart[n]<=gauge_ind[n]<self.nend[n] for n in range(self.ndim)):
                #Set indices relative to this processor's part of grid
                gauge_ind = [gauge_ind[n] - self.nstart[n] for n in range(self.ndim)]
                gauge_path = self.gauge_path+'gauge'+'_'.join(str(coord) for coord in gauge)+'.txt'
                # Is this a good idea???  It should depend on self.overwrite.
                if os.path.isfile(gauge_path): os.remove(gauge_path)
                self.gauges.append(list(gauge_ind))
                self.gauge_files.append(open(gauge_path,'a'))


