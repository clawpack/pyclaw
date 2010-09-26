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

import os
import copy
import logging

import numpy as np

from data import Data
import io

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
    
# Boundary condition user defined function default
def default_user_bc_lower(grid,dim,qbc):
    r"""
    Fills the values of qbc with the correct boundary values
    
    This is a stub function which will return an exception if called.  If you
    want to use a user defined boundary condition replace this function with
    one of your own.
    """
    raise NotImplementedError("Lower user defined boundary condition unimplemented")

def default_user_bc_upper(grid,dim,qbc):
    r"""
    Fills the values of qbc with the correct boundary values
    
    This is a stub function which will return an exception if called.  If you
    want to use a user defined boundary condition replace this function with
    one of your own.
    """
    raise NotImplementedError("Lower user defined boundary condition unimplemented")


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
     - *mthbc_lower* - (int) Lower boundary condition method to be used
     - *mthbc_upper* - (int) Upper boundary condition method to be used
     - *user_bc_lower* - (func) User defined lower boundary condition
     - *user_bc_upper* - (func) User defined upper boundary condition
        
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
    def mthbc_lower():
        doc = r"""(int) - Lower boundary condition fill method. 
                  ``default = 1``"""
        def fget(self): return self._mthbc_lower
        def fset(self,value):
            if value == 2:
                self._mthbc_upper = value
            self._mthbc_lower = value
        return locals()
    mthbc_lower = property(**mthbc_lower())
    _mthbc_lower = 1
    def mthbc_upper():
        doc = r"""(int) - Upper boundary condition fill method. 
                  ``default = 1``"""
        def fget(self): return self._mthbc_upper
        def fset(self,value):
            if value == 2:
                self._mthbc_lower = value
            self._mthbc_upper = value
        return locals()
    mthbc_upper = property(**mthbc_upper())
    _mthbc_upper = 1
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
        self.user_bc_lower = default_user_bc_lower
        r"""(func) - User defined boundary condition function, lower.  
        ``default = None``
        """
        self.user_bc_upper = default_user_bc_upper
        r"""(func) - User defined boundary condition function, upper. 
        ``default = None``"""
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
            
        # Function attribute assignments
    

    def __str__(self):
        output = "Dimension %s" % self.name
        if self.units:
            output += " (%s)" % self.units
        output += ":  (n,d,[lower,upper]) = (%s,%s,[%s,%s])" \
            % (self.n,self.d,self.lower,self.upper)
        output += "  mthbc = (%s,%s)" % (self.mthbc_lower,self.mthbc_upper)
        return output
        

    def qbc_lower(self,grid,qbc):
        r"""
        Apply lower boundary conditions to qbc
        
        Sets the lower coordinate's ghost cells of *qbc* depending on what 
        :attr:`mthbc_lower` is.  If :attr:`mthbc_lower` = 0 then the user 
        boundary condition specified by :attr:`user_bc_lower` is used.  Note 
        that in this case the function :attr:`user_bc_lower` belongs only to 
        this dimension but :attr:`user_bc_lower` could set all user boundary 
        conditions at once with the appropriate calling sequence.
        
        :Input:
         - *grid* - (:class:`Grid`) Grid that the dimension belongs to
         
        :Input/Ouput:
         - *qbc* - (ndarray(...,meqn)) Array with added ghost cells which will
           be set in this routines
        """
        
        # User defined functions
        if self.mthbc_lower == 0:
            self.user_bc_lower(grid,self,qbc)
        # Zero-order extrapolation
        elif self.mthbc_lower == 1:
            qbc[:grid.mbc,...] = qbc[grid.mbc,...]
        # Periodic
        elif self.mthbc_lower == 2:
            qbc[:grid.mbc,...] = qbc[-2*grid.mbc:-grid.mbc,...]
        # Solid wall bc
        elif self.mthbc_lower == 3:
            raise NotImplementedError("Solid wall upper boundary condition not implemented.")
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % self.mthbc_lower)
            
    def qbc_upper(self,grid,qbc):
        r"""
        Apply upper boundary conditions to qbc
        
        Sets the upper coordinate's ghost cells of *qbc* depending on what 
        :attr:`mthbc_upper` is.  If :attr:`mthbc_upper` = 0 then the user 
        boundary condition specified by :attr:`user_bc_upper` is used.  Note 
        that in this case the function :attr:`user_bc_upper` belongs only to 
        this dimension but :attr:`user_bc_upper` could set all user boundary 
        conditions at once with the appropriate calling sequence.
        
        :Input:
         - *grid* - (:class:`Grid`) Grid that the dimension belongs to
         
        :Input/Ouput:
         - *qbc* - (ndarray(...,meqn)) Array with added ghost cells which will
           be set in this routines
        """
        
        # User defined functions
        if self.mthbc_upper == 0:
            self.user_bc_upper(grid,self,qbc)
        # Zero-order extrapolation
        elif self.mthbc_upper == 1:
            qbc[-grid.mbc:,...] = qbc[-grid.mbc-1,...]
        # Periodic
        elif self.mthbc_upper == 2:
            qbc[-grid.mbc:,...] = qbc[grid.mbc:2*grid.mbc,...]
        # Solid wall bc
        elif self.mthbc_upper == 3:
            raise NotImplementedError("Solid wall upper boundary condition not implemented.")
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % self.mthbc_upper)



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
            if self.aux is not None:
                return self.aux.shape[-1]
            return 0
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
    def mthbc_lower():
        doc = r"""(list) - List of lower bc methods"""
        def fget(self): return self.get_dim_attribute('mthbc_lower')
        def fset(self, value): 
            [setattr(getattr(self,dim),'mthbc_lower',value) for dim in self._dimensions]
        return locals()
    def mthbc_upper():
        doc = r"""(list) - List of upper bc methods"""
        def fget(self): return self.get_dim_attribute('mthbc_upper')
        def fset(self, value):
            [setattr(getattr(self,dim),'mthbc_upper',value) for dim in self._dimensions]
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
    mthbc_lower = property(**mthbc_lower())
    mthbc_upper = property(**mthbc_upper())
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
        for dim in self.dimensions:
            if dim.mthbc_lower == 0:
                try:
                    dim.user_bc_lower(self,dim,None)
                except NotImplementedError:
                    logger.debug('Lower BC function for %s has not been set.' % dim.name)
                    valid = False
                except:
                    pass
            if dim.mthbc_upper == 0:
                try:
                    dim.user_bc_upper(self,dim,None)
                except NotImplementedError:
                    logger.debug('Upper BC function for %s has not been set.' % dim.name)
                    valid = False
                except:
                    pass
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
    def empty_q(self,order='C'):
        r"""
        Initialize q to empty
        
        :Input:
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'C'``
        """
        shape = self.get_dim_attribute('n')
        shape.append(self.meqn)
        self.q = np.empty(shape,'d',order=order)
    
    def ones_q(self,order='C'):
        r"""
        Initialize q to all ones
        
        :Input:
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'C'``
        """
        shape = self.get_dim_attribute('n')
        shape.append(self.meqn)
        self.q = np.ones(shape,'d',order=order)
        
    def zeros_q(self,order='C'):
        r"""
        Initialize q to all zeros
        
        :Input:
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'C'``
        """
        shape = self.get_dim_attribute('n')
        shape.append(self.meqn)
        self.q = np.zeros(shape,'d',order=order)
    
    def empty_aux(self,maux,shape=None,order='C'):
        r"""
        Initialize aux to empty with given shape
        
        :Input:
         - *shape* - (tuple) If given, the resulting shape of the auxiliary
           array will be ``shape.append(maux)``.  Otherwise it will be
           ``(dim.n, maux)``
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'C'``
        """
        if shape is None:
            shape = self.n
        shape.append(maux)
        self.aux = np.empty(shape,'d',order=order)
        
    def ones_aux(self,maux,shape=None,order='C'):
        r"""
        Initialize aux to ones with shape
        
        :Input:
         - *shape* - (tuple) If given, the resulting shape of the auxiliary
           array will be ``shape.append(maux)``.  Otherwise it will be
           ``(dim.n, maux)``
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'C'``
        """
        if shape is None:
            shape = self.n
        shape.append(maux)
        self.aux = np.ones(shape,'d',order=order)
        
    def zeros_aux(self,maux,shape=None,order='C'):
        r"""
        Initialize aux to zeros with shape
        
        :Input:
         - *shape* - (tuple) If given, the resulting shape of the auxiliary
           array will be ``shape.append(maux)``.  Otherwise it will be
           ``(dim.n, maux)``
         - *order* - (string) Order of array, must be understood by numpy
           ``default = 'C'``
        """
        if shape is None:
            shape = self.n
        shape.append(maux)
        self.aux = np.zeros(shape,'d',order=order)
        
                
    def qbc(self):
        r"""
        Appends boundary conditions to q
    
        This function returns an array of dimension determined by the 
        :attr:`mbc` attribute.  The type of boundary condition set is 
        determined by :attr:`mthbc_lower` and :attr:`mthbc_upper` for the 
        approprate dimension.  Valid values for :attr:`mthbc_lower` and 
        :attr:`mthbc_upper` include:
        
        0. A user defined boundary condition will be used, the appropriate 
           Dimension method user_bc_lower or user_bc_upper will be called
        1. Zero-order extrapolation
        2. Periodic boundary conditions
        3. Wall boundary conditions, it is assumed that the second equation
           represents velocity or momentum
    
        :Output:
         - (ndarray(...,meqn)) q array with boundary ghost cells added and set
         
        .. note:: 
            Note that for user defined boundary conditions, the array sent to
            the boundary condition has not been rolled. 
        """
        
        # Create ghost cell array
        dim_string = ','.join( ('2*self.mbc+self.%s.n' % name for name in self._dimensions) )
        exec("qbc = np.zeros( (%s,self.meqn) )" % dim_string)
        dim_string = ','.join( ('self.mbc:-self.mbc' for name in self._dimensions) )
        exec("qbc[%s,:] = self.q" % dim_string)
        
        for i in xrange(len(self._dimensions)):
            dim = getattr(self,self._dimensions[i])
            
            # If a user defined boundary condition is being used, send it on,
            # otherwise roll the axis to front position and operate on it
            if dim.mthbc_lower == 0:
                dim.qbc_lower(self,qbc)
            else:
                dim.qbc_lower(self,np.rollaxis(qbc,i))
            if dim.mthbc_upper == 0:
                dim.qbc_upper(self,qbc)
            else:
                dim.qbc_upper(self,np.rollaxis(qbc,i))
            
        return qbc


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



# ============================================================================
#  Solution Class
# ============================================================================
class Solution(object):
    r"""
    Pyclaw grid container class
        
    :Input and Output:
    
        Input and output of solution objects is handle via the io package.
        Solution contains the generic methods :meth:`write`, :meth:`read` and
        :meth:`plot` which then figure out the correct method to call.  Please
        see the io package for the particulars of each format and method and 
        the methods in this class for general input and output information.
    
    :Properties:
    
        If there is only one grid belonging to this grid, the solution will
        appear to have many of the attributes assigned to its one grid.  Some
        parameters that have in the past been parameters for all grids are
        also reachable although Solution does not check to see if these
        parameters are truly universal.

        Grid Attributes:
            't','meqn','mbc','q','aux','capa','aux_global','dimensions'
            
    :Initialization:
        
        The initialization of a Solution can happen on of these ways
            1. args is empty and an empty Solution is created
            2. args is a single Grid or list of Grids
            3. args is a single Dimension or list of Dimensions
            4. args is a variable number of arguments that describes the 
               location of a file to be read in to initialize the object
            5. args is a data object with the corresponding data fields
        
        Input:
            - if args == () -> Empty Solution object
            - if args == Grids -> Grids are appended to grids list
            - if args == Dimensions -> A single Grid with the given
              Dimensions is created and appended to the grids list
            - if args == frame, format='ascii',path='./',file_prefix='fort'
            - if args == Data, Create a new single grid solution based off of 
              what is in args.
    
    """

    # ========== Attributes ==================================================
    
    # ========== Properties ==================================================
    def grid():
        doc = r"""(:class:`Grid`) - Base grid is returned"""
        def fget(self): return self.grids[0]
        return locals()
    def t():
        doc = r"""(float) - :attr:`Grid.t` of base grid"""
        def fget(self): return self._get_base_grid_attribute('t')
        def fset(self, value): self.set_all_grids('t',value)
        return locals()
    def meqn():
        doc = r"""(int) - :attr:`Grid.meqn` of base grid"""
        def fget(self): return self._get_base_grid_attribute('meqn')
        def fset(self, value): self.set_all_grids('meqn',value)
        return locals()   
    def mbc():
        doc = r"""(int) - :attr:`Grid.mbc` of base grid"""
        def fget(self): return self._get_base_grid_attribute('mbc')
        def fset(self, value): self.set_all_grids('mbc',value)
        return locals()   
    def q():
        doc = r"""(ndarray(...,:attr:`Grid.meqn`)) - :attr:`Grid.q` of base 
                  grid"""
        def fget(self): return self._get_base_grid_attribute('q')
        return locals()
    def aux():
        doc = r"""(ndarray(...,:attr:`Grid.maux`)) - :attr:`Grid.aux` of base 
                  grid"""
        def fget(self): return self._get_base_grid_attribute('aux')
        def fset(self, value): 
            if len(self.grids) == 1: 
                setattr(self.grids[0],'aux',value)
        return locals()  
    def capa():
        doc = r"""(ndarray(...)) - :attr:`Grid.capa` of base grid"""
        def fget(self): return self._get_base_grid_attribute('capa')
        def fset(self, value):
            if len(self.grids) == 1: 
                setattr(self.grids[0],'capa',value)
        return locals()  
    def aux_global():
        doc = r"""(dict) - :attr:`Grid.aux_global` of base grid"""
        def fget(self): return self._get_base_grid_attribute('aux_global')
        def fset(self, value):
            if len(self.grids) == 1: 
                setattr(self.grids[0],'aux_global',value)
        return locals()
    def maux():
        doc = r"""(int) - :attr:`Grid.maux` of base grid"""
        def fget(self): return self._get_base_grid_attribute('maux')
        return locals()
    def ndim():
        doc = r"""(int) - :attr:`Grid.ndim` of base grid"""
        def fget(self): return self._get_base_grid_attribute('ndim')
        return locals()
    def dimensions():
        doc = r"""(list) - :attr:`Grid.dimensions` of base grid"""
        def fget(self): return self._get_base_grid_attribute('dimensions')
        return locals()
    def n():
        doc = r"""(list) - :attr:`Grid.n` of base grid"""
        def fget(self): return self._get_base_grid_attribute('n')
        return locals()
    def name():
        doc = r"""(list) - :attr:`Grid.name` of base grid"""
        def fget(self): return self._get_base_grid_attribute('name')
        return locals()
    def lower():
        doc = r"""(list) - :attr:`Grid.lower` of base grid"""
        def fget(self): return self._get_base_grid_attribute('lower')
        return locals()
    def upper():
        doc = r"""(list) - :attr:`Grid.upper` of base grid"""
        def fget(self): return self._get_base_grid_attribute('upper')
        return locals()
    def d():
        doc = r"""(list) - :attr:`Grid.d` of base grid"""
        def fget(self): return self._get_base_grid_attribute('d')
        return locals()
    def units():
        doc = r"""(list) - :attr:`Grid.units` of base grid"""
        def fget(self): return self._get_base_grid_attribute('units')
        return locals()
    def mthbc_lower():
        doc = r"""(list) - :attr:`Grid.mthbc_lower` of base grid"""
        def fget(self): return self._get_base_grid_attribute('mthbc_lower')
        def fset(self, value): self.set_all_grids('mthbc_lower',value)
        return locals()
    def mthbc_upper():
        doc = r"""(list) - :attr:`Grid.mthbc_upper` of base grid"""
        def fget(self): return self._get_base_grid_attribute('mthbc_upper')
        def fset(self, value): self.set_all_grids('mthbc_upper',value)
        return locals()
    def center():
        doc = r"""(list) - :attr:`Grid.center` of base grid"""
        def fget(self): return self._get_base_grid_attribute('center')
        return locals()
    def edge():
        doc = r"""(list) - :attr:`Grid.edge` of base grid"""
        def fget(self): return self._get_base_grid_attribute('edge')
        return locals()
    def p_center():
        doc = r"""(list) - :attr:`Grid.p_center` of base grid"""
        def fget(self): return self._get_base_grid_attribute('p_center')
        return locals()
    def p_edge():
        doc = r"""(list) - :attr:`Grid.p_edge` of base grid"""
        def fget(self): return self._get_base_grid_attribute('p_edge')
        return locals()
    def c_center():
        doc = r"""(list) - :attr:`Grid.c_center` of base grid"""
        def fget(self): return self._get_base_grid_attribute('c_center')
        return locals()
    def c_edge():
        doc = r"""(list) - :attr:`Grid.c_edge` of base grid"""
        def fget(self): return self._get_base_grid_attribute('c_edge')
        return locals()
        
    grid = property(**grid())
    t = property(**t())
    meqn = property(**meqn()) 
    mbc = property(**mbc())
    q = property(**q())
    aux = property(**aux())
    capa = property(**capa())
    aux_global = property(**aux_global())
    maux = property(**maux())
    ndim = property(**ndim())
    dimensions = property(**dimensions())
    n = property(**n())
    name = property(**name())
    lower = property(**lower())
    upper = property(**upper())
    d = property(**d())
    units = property(**units())
    mthbc_lower = property(**mthbc_lower())
    mthbc_upper = property(**mthbc_upper())
    center = property(**center())
    edge = property(**edge())
    p_center = property(**p_center())
    p_edge = property(**p_edge())
    c_center = property(**c_center())
    c_edge = property(**c_edge())
    

    # ========== Class Methods ===============================================
    def __init__(self,*arg,**kargs):
        r"""Solution Initiatlization Routine
        
        See :class:`Solution` for more info.
        """
        self.grids = []
        r"""(list) - List of grids in this solution"""
        
        # If arg is non-zero, we are reading in a solution, otherwise, we
        # create an empty Solution
        if len(arg) > 0:
            # Single Grid
            if isinstance(arg[0],Grid):
                self.grids.append(arg[0])
            # Single Dimension
            elif isinstance(arg[0],Dimension):
                self.grids.append(Grid(arg[0]))
            elif isinstance(arg[0],list):
                # List of Grids
                if isinstance(arg[0][0],Grid):
                    self.grids = arg[0]
                # List of Dimensions
                elif isinstance(arg[0][0],Dimension):
                    self.grids.append(Grid(arg[0]))
                else:
                    raise Exception("Invalid argument list")
            elif isinstance(arg[0],int): 
                frame = arg[0]
                defaults = {'format':'ascii','path':'./','file_prefix':None,
                    'read_aux':True,'options':{}}
   
                for (k,v) in defaults.iteritems():    
                    if kargs.has_key(k):
                        exec("%s = kargs['%s']" % (k,k))
                    else:
                        exec('%s = v' % k)
                self.read(frame,path,format,file_prefix,read_aux,options)
            elif isinstance(arg[0],Data):
                data = arg[0] 
                # Create dimensions
                if data.ndim == 1:
                    x = Dimension('x',data.xlower,data.xupper,data.mx,
                                    mthbc_lower=data.mthbc_xlower,
                                    mthbc_upper=data.mthbc_xupper)
                    grid = Grid([x])
                elif data.ndim == 2:
                    x = Dimension('x',data.xlower,data.xupper,data.mx,
                                    mthbc_lower=data.mthbc_xlower,
                                    mthbc_upper=data.mthbc_xupper)
                    y = Dimension('y',data.ylower,data.yupper,data.my,
                                    mthbc_lower=data.mthbc_ylower,
                                    mthbc_upper=data.mthbc_yupper)
                    grid = Grid([x,y])
                elif data.ndim == 3:
                    x = Dimension('x',data.xlower,data.xupper,data.mx,
                                    mthbc_lower=data.mthbc_xlower,
                                    mthbc_upper=data.mthbc_xupper)
                    y = Dimension('y',data.ylower,data.yupper,data.my,
                                    mthbc_lower=data.mthbc_ylower,
                                    mthbc_upper=data.mthbc_yupper)
                    z = Dimension('z',data.zlower,data.zupper,data.mz,
                                    mthbc_lower=data.mthbc_zlower,
                                    mthbc_upper=data.mthbc_zupper)
                    grid = Grid([x,y,z])
                
                # General grid properties
                grid.mbc = data.mbc
                grid.t = data.t0
                grid.meqn = data.meqn
                
                # Add grid to solution
                self.grids.append(grid)
            else:
                raise Exception("Invalid argument list")
                
                
    def is_valid(self):
        r"""
        Checks to see if this solution is valid
        
        The Solution checks to make sure it is valid by checking each of its
        grids.  If an invalid grid is found, a message is logged what
        specifically made this solution invalid.
        
        :Output:
         - (bool) - True if valid, false otherwise
        """
        return reduce(lambda x,y: x*y,[grid.is_valid() for grid in self.grids])


    def __str__(self):
        output = "Grids:\n"
        # This is information about each of the grids
        for grid in self.grids:
            output = ''.join((output,'  %s:\n' % grid.gridno))
            output = output + str(grid)
        return str(output)
    
    
    def set_all_grids(self,attr,value,overwrite=True):
        r"""
        Sets all member grids attribute 'attr' to value
        
        :Input:
         - *attr* - (string) Attribute name to be set
         - *value* - (id) Value for attribute
         - *overwrite* - (bool) Whether to overwrite the attribute if it 
           already exists.  ``default = True``
        """
        [setattr(grid,attr,value) for grid in self.grids if getattr(grid,attr)
                    is None or overwrite]
    
                    
    def _get_base_grid_attribute(self, name):
        r"""
        Return base grid attribute name
        
        :Output:
         - (id) - Value of attribute from ``grids[0]``
        """
        return getattr(self.grids[0],name)
    
    
    def __copy__(self):
        return self.__class__(self)
    
    
    def __deepcopy__(self, memo={}):
        # Create basic container
        result = self.__class__()
        result.__init__()
        
        # Populate the grids
        for grid in self.grids:
            result.grids.append(copy.deepcopy(grid))
        
        return result
    
    
    # ========== IO Functions ================================================
    def write(self,frame,path='./',format='ascii',file_prefix=None,
                write_aux=False,options={}):
        r"""
        Write out a representation of the solution

        Writes out a suitable representation of this solution object based on
        the format requested.  The path is built from the optional path and
        file_prefix arguments.  Will raise an IOError if unsuccessful.

        :Input:
         - *frame* - (int) Frame number to append to the file output
         - *path* - (string) Root path, will try and create the path if it 
           does not already exist. ``default = './'``
         - *format* - (string or list of strings) a string or list of strings 
           containing the desired output formats. ``default = 'ascii'``
         - *file_prefix* - (string) Prefix for the file name.  Defaults to
           the particular io modules default.
         - *write_aux* - (book) Write the auxillary array out as well if 
           present. ``default = False``
         - *options* - (dict) Dictionary of optional arguments dependent on 
           which format is being used. ``default = {}``
        """
        # Determine if we need to create the path
        path = os.path.expandvars(os.path.expanduser(path))
        if not os.path.exists(path):
            # Attempt to constuct the path  
            os.makedirs(path)

        # Call the correct write function based on the output format
        if isinstance(format,str):
            format_list = [format]
        elif isinstance(format,list):
            format_list = format
        # Loop over list of formats requested
        for form in format_list:
            write_func = eval('io.write_%s' % form)
            if file_prefix is None:
                write_func(self,frame,path,write_aux=write_aux,
                            options=options)
            else:
                write_func(self,frame,path,file_prefix=file_prefix,
                                write_aux=write_aux,options=options)
            msg = "Wrote out solution in format %s for time t=%s" % (form,self.t)
            logging.getLogger('io').info(msg)
        
        
    def read(self,frame,path='./',format='ascii',file_prefix=None,
                read_aux=False,options={}):
        r"""
        Reads in a Solution object from a file
        
        Reads in and initializes this Solution with the data specified.  This 
        function will raise an IOError if it was unsuccessful.  

        Any format must conform to the following call signiture and return
        True if the file has been successfully read into the given solution or
        False otherwise.  Options is a dictionary of parameters that each
        format can specify.  See the ascii module for an example.::
        
            read_<format>(solution,path,frame,file_prefix,options={})
            
        ``<format>`` is the name of the format in question.
        
        :Input:
         - *frame* - (int) Frame number to be read in
         - *path* - (string) Base path to the files to be read. 
           ``default = './'``
         - *format* - (string) Format of the file, should match on of the 
           modules inside of the io package.  ``default = 'ascii'``
         - *file_prefix* - (string) Name prefix in front of all the files, 
           defaults to whatever the format defaults to, e.g. fort for ascii
         - *options* - (dict) Dictionary of optional arguments dependent on 
           the format being read in.  ``default = {}``
            
        :Output:
         - (bool) - True if read was successful, False otherwise
        """
        
        path = os.path.expandvars(os.path.expanduser(path))
        read_func = eval('io.read_%s' % format)
        if file_prefix is None:
            read_func(self,frame,path,read_aux=read_aux,options=options)
        else:
            read_func(self,frame,path,file_prefix=file_prefix,
                                    read_aux=read_aux,options=options)
        logging.getLogger('io').info("Read in solution for time t=%s" % self.t)
        
        
    def plot(self):
        r"""
        Plot the solution
        """
        raise NotImplementedError("Direct solution plotting has not been " +
            "implemented as of yet, please refer to the plotting module for" +
            " how to plot solutions.")
