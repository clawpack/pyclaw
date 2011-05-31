#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing all Pyclaw solution objects

:Authors:
    Kyle T. Mandli (2008-08-07) Initial version
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

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
                    self.edge[i] = self.lower + i*self.d 
            return self._edge
        return locals()
    edge = property(**edge())
    _edge = None
    def centerghost():
        doc = r"""(ndarrary(:)) - Location of all grid cell center coordinates
        for this dimension, including ghost cells"""
        def fget(self): 
            if self._centerghost is None:
                self._centerghost = np.empty(self.n+2*self.mbc)
                for i in xrange(0-self.mbc,self.n+self.mbc):
                    self.centerghost[i] = self.lower + (i+0.5)*self.d
            return self._centerghost
        return locals()
    centerghost = property(**centerghost())
    _centerghost = None
    def center():
        doc = r"""(ndarrary(:)) - Location of all grid cell center coordinates
        for this dimension"""
        def fget(self): 
            if self._center is None:
                self._center = np.empty(self.n)
                for i in xrange(0,self.n):
                    self.center[i] = self.lower + (i+0.5)*self.d
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
        self.n = 1
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
    
        Each grid has a value for :attr:`level`, :attr:`gridno`, :attr:`t`, 
        :attr:`mbc`, :attr:`meqn` and :attr:`aux_global`.  These correspond
        to global grid traits and determine many of the properties and sizes
        of the data arrays.
        
    :Grid Data:
    
        The arrays :attr:`q`, :attr:`aux` and :attr:`capa` have variable 
        extents based on the set of dimensions present and the values of 
        :attr:`meqn` and :attr:`maux`.  Note that these are initialy set to 
        None so need to be instantiated.  For convenience, the methods
        :meth:`emtpy_q`, :meth:`ones_q`, and :meth:`zeros_q` for ``q`` and
        :meth:`emtpy_aux`, :meth:`ones_aux`, and :meth:`zeros_aux` for ``aux``
        are provided to initialize these arrays.  The :attr:`capa` array is 
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
    def maux():
        doc = r"""(int) - Rank of auxiliary array"""
        def fget(self):
            if self.aux is not None: return self.aux.shape[0]
            else: return 0
        return locals()
    def n():
        doc = r"""(list) - List of the number of grid cells in each dimension"""
        def fget(self): return self.get_dim_attribute('n')
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
    maux = property(**maux())
    n = property(**n())
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
        self.t = 0.0
        r"""(float) - Current time represented on this grid, 
            ``default = 0.0``"""
        self.mbc = 2
        r"""(int) - Number of ghost cells along the boundaries, 
            ``default = 2``"""
        self.meqn = 1
        r"""(int) - Dimension of q array for this grid, ``default = 1``"""
        self.q = None
        r"""(ndarray(...,meqn)) - Cell averaged quantity being evolved."""
        self.aux = None
        r"""(ndarray(...,maux)) - Auxiliary array for this grid containing per 
            cell information"""
        self.capa = None
        r"""(ndarray(...)) - Capacity array for this grid, ``default = 1.0``"""
        self.aux_global = {}
        r"""(dict) - Dictionary of global values for this grid, 
            ``default = {}``"""
        self.mapc2p = default_mapc2p
        r"""(func) - Grid mapping function"""
        
        # Dimension parsing
        if isinstance(dimensions,Dimension):
            dimensions = [dimensions]
        self._dimensions = []
        for dim in dimensions:
            self.add_dimension(dim)
        for dim in self.dimensions:
            dim.mbc=self.mbc
    
    
    def __str__(self):
        output = "Grid %s:\n" % self.gridno
        output += "  t=%s mbc=%s meqn=%s\n  " % (self.t,self.mbc,self.meqn)
        output += '\n  '.join((str(getattr(self,dim)) for dim in self._dimensions))
        output += '\n'
        if self.q is not None:
            output += "  q.shape=%s" % str(self.q.shape)
        if self.aux is not None:
            output += " aux.shape=%s" % str(self.aux.shape)
        if self.capa is not None:
            output += " capa.shape=%s" % str(self.capa.shape)
        return output
    
    
    def is_valid(self):
        r"""
        Checks to see if this grid is valid
        
        The grid is declared valid based on the following criteria:
            - :attr:`q` is not None
            - All dimensions with boundary conditions set to 0 have valid bc 
              functions
            
        A debug logger message will be sent documenting exactly what was not 
        valid.
            
        :Output:
         - (bool) - True if valid, false otherwise.
        
        """
        valid = True
        logger = logging.getLogger('solution')
        if self.q is None:
            logger.debug('The array q has not been initialized.')
            valid = False
        return valid
    
    
    def set_aux_global(self,data_path):
        r"""
        Convenience routine for parsing data files into the aux_global dict
        
        Puts the data from the file at path and puts it into the aux_global
        dictionary using :class:`Data`.  This assumes that all values in the
        file are to be put into the :attr:`aux_global` dictionary and the file
        conforms to that which :class:`Data` can read.
        
        :Input:
         - *data_path* - (string) Path to the file to be read in
         
        .. note:: 
            This will not replace the :attr:`aux_global` dictionary but add to
            or replace any values already in it.
        """
        data = Data(data_files=data_path)
        for (k,v) in data.iteritems():
            self.aux_global[k] = v
        
    def set_cparam(self,fortran_module):
        """
        Set the variables in fortran_module.cparam to the corresponding values in
        grid.aux_global.  This is the mechanism for passing scalar variables to the
        Fortran Riemann solvers; cparam must be defined as a common block in the
        Riemann solver.

        This function should be called from solver.setup().  This seems like a fragile
        interdependency between solver and grid; perhaps aux_global should belong
        to solver instead of grid.

        This function also checks that the set of variables defined in cparam 
        all appear in aux_global.

        """
        if hasattr(fortran_module,'cparam'):
            if not set(dir(fortran_module.cparam)) <= set(self.aux_global.keys()):
                raise Exception('Some required value(s) in the cparam common block in the Riemann solver have not been set in aux_global.')
            for global_var_name,global_var_value in self.aux_global.iteritems(): 
                setattr(fortran_module.cparam,global_var_name,global_var_value)

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
        
        
    def remove_dimension(self, name):
        r"""
        Remove the dimension named name
        
        :Input:
         - *name* - (string) removes the dimension named name
        """
        
        # Remove coordinate array
        self._p_center.pop(self._dimension.index(name))
        self._c_center.pop(self._dimension.index(name))
        self._p_edge.pop(self._dimension.index(name))
        self._c_edge.pop(self._dimension.index(name))
        
        # Delete the dimension
        self._dimensions.remove(name)
        exec('del(self.%s)' % name)
    
    
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
        
        for attr in ('level','gridno','t','mbc','meqn','_p_center','_p_edge',
                        '_c_center','_c_edge'):
            setattr(result,attr,copy.deepcopy(getattr(self,attr)))
        
        if self.q is not None:
            result.q = copy.deepcopy(self.q)
        if self.aux is not None:
            result.aux = copy.deepcopy(self.aux)
        if self.capa is not None:
            result.capa = copy.deepcopy(self.capa)
        result.aux_global = copy.deepcopy(self.aux_global)
        
        result.mapc2p = self.mapc2p
        
        return result
        
    
    # ========== Grid Operations =============================================
    # Convenience routines for initialization of q and aux
    def empty_q(self,order='F'):
        r"""
        Initialize q to empty
        
        :Input:
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'F'``
        """
        shape = [dim.nend-dim.nstart for dim in self.dimensions]
        shape.insert(0,self.meqn)
        self.q = np.empty(shape,'d',order=order)
    
    def ones_q(self,order='F'):
        r"""
        Initialize q to all ones
        
        :Input:
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'F'``
        """
        shape = [dim.nend-dim.nstart for dim in self.dimensions]
        shape.insert(0,self.meqn)
        self.q = np.ones(shape,'d',order=order)
        
    def zeros_q(self,order='F'):
        r"""
        Initialize q to all zeros
        
        :Input:
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'F'``
        """
        shape = [dim.nend-dim.nstart for dim in self.dimensions]
        shape.insert(0,self.meqn)
        self.q = np.zeros(shape,'d',order=order)
    
    def empty_aux(self,maux,shape=None,order='F'):
        r"""
        Initialize aux to empty with given shape
        
        :Input:
         - *shape* - (tuple) If given, the resulting shape of the auxiliary
           array will be ``shape.append(maux)``.  Otherwise it will be
           ``(dim.n, maux)``
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'F'``
        """
        if shape is None:
            shape = self.n
        shape.insert(0,maux)
        self.aux = np.empty(shape,'d',order=order)
        
    def ones_aux(self,maux,shape=None,order='F'):
        r"""
        Initialize aux to ones with shape
        
        :Input:
         - *shape* - (tuple) If given, the resulting shape of the auxiliary
           array will be ``shape.append(maux)``.  Otherwise it will be
           ``(dim.n, maux)``
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'F'``
        """
        if shape is None:
            shape = self.n
        shape.insert(0,maux)
        self.aux = np.ones(shape,'d',order=order)
        
    def zeros_aux(self,maux,shape=None,order='F'):
        r"""
        Initialize aux to zeros with shape
        
        :Input:
         - *shape* - (tuple) If given, the resulting shape of the auxiliary
           array will be ``shape.append(maux)``.  Otherwise it will be
           ``(dim.n, maux)``
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'F'``
        """
        if shape is None:
            shape = self.n
        shape.insert(0,maux)
        self.aux = np.zeros(shape,'d',order=order)

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
            self._p_center = [None for i in xrange(self.ndim)]

            # Special case
            if self.ndim == 1:
                self._p_center[0] = self.mapc2p(self,self.dimensions[0].center)
            # Higer dimensional calculate center arrays
            else:    
                # Produce ndim mesh grid function
                mgrid = np.lib.index_tricks.nd_grid()
            
                # Calculate length of full arrays and indexing string
                array_length = reduce(lambda x,y: x*y, self.get_dim_attribute('n'))
                index_str = ','.join(('0:self.%s.n' % dim for dim in self._dimensions))

                # Produce index arrays
                exec('index = mgrid[%s]' % index_str)
                index_str = ','.join((':' for dim in self._dimensions))
            
                # Create meshgrid like arrays
                array_list = []
                i = 0
                for center_array in self.get_dim_attribute('center'):
                    exec("array_list.append(np.reshape(center_array[index[%i,%s]],(array_length)))" % (i,index_str))
                    i += 1
            
                # Actual mapping
                p = self.mapc2p(self,array_list)
            
                # Reshape arrays for storage
                for i in xrange(self.ndim):
                    self._p_center[i] = np.reshape(p[i], 
                                                self.get_dim_attribute('n'))


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
                # Produce ndim mesh grid function
                mgrid = np.lib.index_tricks.nd_grid()
            
                # Calculate length of full arrays and indexing string
                array_length = reduce(lambda x,y: x*y, [n+1 for n in self.get_dim_attribute('n')])
                index_str = ','.join(('0:self.%s.n+1' % dim for dim in self._dimensions))

                # Produce index arrays
                exec('index = mgrid[%s]' % index_str)
                index_str = ','.join((':' for dim in self._dimensions))
            
                # Create meshgrid like arrays
                array_list = []
                i = 0
                for edge_array in self.get_dim_attribute('edge'):
                    exec("array_list.append(np.reshape(edge_array[index[%i,%s]],(array_length)))" % (i,index_str))
                    i += 1
            
                # Actual mapping
                p = self.mapc2p(self,array_list)
            
                # Reshape arrays for storage
                for i in xrange(self.ndim):
                    self._p_edge[i] = np.reshape(p[i], 
                                [n+1 for n in self.get_dim_attribute('n')])
                                

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
            self._c_center = [None for i in xrange(self.ndim)]
            
            # For one dimension, the center and edge arrays are equivalent
            if self.ndim == 1:
                self._c_center[0] = self.dimensions[0].center
            else:
                # Produce ndim mesh grid function
                mgrid = np.lib.index_tricks.nd_grid()
            
                # Create index arrays
                index_str = ','.join(('0:self.%s.n' % dim for dim in self._dimensions))
                exec('index = mgrid[%s]' % index_str)
            
                # Create c_center array
                index_str = ','.join((':' for dim in self._dimensions))
                for i in xrange(self.ndim):
                    exec("self._c_center[i] = self.dimensions[i].center[index[i,%s]]" % index_str)
                    
                    
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
                # Produce ndim mesh grid function
                mgrid = np.lib.index_tricks.nd_grid()
            
                # Create index arrays
                index_str = ','.join(('0:self.%s.n+1' % dim for dim in self._dimensions))
                exec('index = mgrid[%s]' % index_str)
            
                # Create c_center array
                index_str = ','.join((':' for dim in self._dimensions))
                for i in xrange(self.ndim):
                    exec("self._c_edge[i] = self.dimensions[i].edge[index[i,%s]]" % index_str)



