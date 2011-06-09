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
from state import State
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
    
        If there is only one state and grid belonging to this solution, 
        the solution will appear to have many of the attributes assigned to its
        one state and grid.  Some parameters that have in the past been
        parameters for all grids are also reachable although Solution does not
        check to see if these parameters are truly universal.

        Grid Attributes:
            'dimensions'
        State Attributes:
            't','meqn','q','aux','capa','aux_global'
            
    :Initialization:
        
        The initialization of a Solution can happen on of these ways
            1. args is empty and an empty Solution is created
            2. args is a single State or list of States
            2. args is a single Grid or list of Grids
            3. args is a single Dimension or list of Dimensions
            4. args is a variable number of arguments that describes the 
               location of a file to be read in to initialize the object
            5. args is a data object with the corresponding data fields
        
        Input:
            - if args == () -> Empty Solution object
            - if args == States -> States are appended to states list
            - if args == Grids -> States are initialized with these Grids 
              and appended to states list
            - if args == Dimensions -> A single Grid with the given
              Dimensions is created, a state is initalized with this Grid
              and appended to the states list
            - if args == frame, format='ascii',path='./',file_prefix='fort'
            - if args == Data, Create a new single grid solution based off of 
              what is in args.
    
    """

    # ========== Attributes ==================================================
    
    # ========== Properties ==================================================
    def state():
        doc = r"""(:class:`State`) - Base state is returned"""
        def fget(self): return self.states[0]
        return locals()
    def grid():
        doc = r"""(:class:`Grid`) - Base state's grid is returned"""
        def fget(self): return self.states[0].grid
        return locals()
    def t():
        doc = r"""(float) - :attr:`State.t` of base state"""
        def fget(self): return self._get_base_state_attribute('t')
        def fset(self, value): self.set_all_states('t',value)
        return locals()
    def meqn():
        doc = r"""(int) - :attr:`State.meqn` of base state"""
        def fget(self): return self._get_base_state_attribute('meqn')
        def fset(self, value): self.set_all_states('meqn',value)
        return locals()   
    def q():
        doc = r"""(ndarray(...,:attr:`State.meqn`)) - :attr:`State.q` of base 
                  state"""
        def fget(self): return self._get_base_state_attribute('q')
        return locals()
    def aux():
        doc = r"""(ndarray(...,:attr:`State.maux`)) - :attr:`State.aux` of base 
                  state"""
        def fget(self): return self._get_base_state_attribute('aux')
        def fset(self, value): 
            if len(self.states) == 1: 
                setattr(self.states[0],'aux',value)
        return locals()  
    def capa():
        doc = r"""(ndarray(...)) - :attr:`State.capa` of base state"""
        def fget(self): return self._get_base_state_attribute('capa')
        def fset(self, value):
            if len(self.states) == 1: 
                setattr(self.states[0],'capa',value)
        return locals()  
    def aux_global():
        doc = r"""(dict) - :attr:`State.aux_global` of base state"""
        def fget(self): return self._get_base_state_attribute('aux_global')
        def fset(self, value):
            if len(self.states) == 1: 
                setattr(self.states[0],'aux_global',value)
        return locals()
    def maux():
        doc = r"""(int) - :attr:`State.maux` of base state"""
        def fget(self): return self._get_base_state_attribute('maux')
        return locals()
    def ndim():
        doc = r"""(int) - :attr:`Grid.ndim` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('ndim')
        return locals()
    def dimensions():
        doc = r"""(list) - :attr:`Grid.dimensions` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('dimensions')
        return locals()
    def n():
        doc = r"""(list) - :attr:`Grid.n` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('n')
        return locals()
    def name():
        doc = r"""(list) - :attr:`Grid.name` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('name')
        return locals()
    def lower():
        doc = r"""(list) - :attr:`Grid.lower` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('lower')
        return locals()
    def upper():
        doc = r"""(list) - :attr:`Grid.upper` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('upper')
        return locals()
    def d():
        doc = r"""(list) - :attr:`Grid.d` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('d')
        return locals()
    def units():
        doc = r"""(list) - :attr:`Grid.units` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('units')
        return locals()
    def center():
        doc = r"""(list) - :attr:`Grid.center` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('center')
        return locals()
    def edge():
        doc = r"""(list) - :attr:`Grid.edge` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('edge')
        return locals()
    def p_center():
        doc = r"""(list) - :attr:`Grid.p_center` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('p_center')
        return locals()
    def p_edge():
        doc = r"""(list) - :attr:`Grid.p_edge` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('p_edge')
        return locals()
    def c_center():
        doc = r"""(list) - :attr:`Grid.c_center` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('c_center')
        return locals()
    def c_edge():
        doc = r"""(list) - :attr:`Grid.c_edge` of base state.grid"""
        def fget(self): return self._get_base_grid_attribute('c_edge')
        return locals()
        
    state = property(**state())
    grid = property(**grid())
    t = property(**t())
    meqn = property(**meqn()) 
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
        self.states = []
        r"""(list) - List of grids in this solution"""
        # If arg is non-zero, we are reading in a solution, otherwise, we
        # create an empty Solution
        if len(arg) > 0:
            # Single State
            if isinstance(arg[0],State):
                self.states.append(arg[0])
            # Single Grid
            elif isinstance(arg[0],Grid):
                self.states.append(State(arg[0]))
            # Single Dimension
            elif isinstance(arg[0],Dimension):
                self.states.append(State(Grid(arg[0])))
            elif isinstance(arg[0],list):
                # List of States
                if isinstance(arg[0][0],State):
                    self.states = arg[0]
                # List of Grids
                if isinstance(arg[0][0],Grid):
                    self.states = Grid(arg[0])
                # List of Dimensions
                elif isinstance(arg[0][0],Dimension):
                    self.state.append(State(Grid(arg[0])))
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
                
                state = State(grid)
                # General grid properties
                state.t = data.t0
                state.meqn = data.meqn
                
                # Add grid to solution
                self.states.append(state)
            else:
                raise Exception("Invalid argument list")
                
                
    def is_valid(self):
        r"""
        Checks to see if this solution is valid
        
        The Solution checks to make sure it is valid by checking each of its
        states.  If an invalid state is found, a message is logged what
        specifically made this solution invalid.
       
        :Output:
         - (bool) - True if valid, false otherwise
        """
        return all([state.is_valid() for state in self.states])


    def __str__(self):
        output = "states:\n"
        # This is information about each of the states
        for state in self.states:
            output = ''.join((output,'  %s:\n' % state.stateno))
            output = output + str(state)
        return str(output)
    
    
    def set_all_states(self,attr,value,overwrite=True):
        r"""
        Sets all member states attribute 'attr' to value
        
        :Input:
         - *attr* - (string) Attribute name to be set
         - *value* - (id) Value for attribute
         - *overwrite* - (bool) Whether to overwrite the attribute if it 
           already exists.  ``default = True``
        """
        [setattr(state,attr,value) for state in self.states if getattr(state,attr)
                    is None or overwrite]
    
                    
    def _get_base_state_attribute(self, name):
        r"""
        Return base state attribute
        
        :Output:
         - (id) - Value of attribute from ``states[0]``
        """
        return getattr(self.states[0],name)
    
    def _get_base_grid_attribute(self, name):
        r"""
        Return base state.grid attribute name
        
        :Output:
         - (id) - Value of attribute from ``states[0].grid``
        """
        return getattr(self.states[0].grid,name)
    
    
    def __copy__(self):
        return self.__class__(self)
    
    
    def __deepcopy__(self, memo={}):
        # Create basic container
        result = self.__class__()
        result.__init__()
        
        # Populate the states
        for state in self.states:
            result.states.append(copy.deepcopy(state))
        
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
