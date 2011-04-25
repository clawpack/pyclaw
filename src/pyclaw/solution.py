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

from data import Data
from grid import Grid, Dimension
import io

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
    
        If there is only one grid belonging to this solution, the solution will
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
            try:
                os.makedirs(path)
            except OSError:
                print "directory already exists, ignoring"  

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
