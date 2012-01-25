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
import logging

from pyclaw.grid import Grid, Dimension
from pyclaw.state import State
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
            't','num_eqn','q','aux','capa','problem_data'
            
    :Initialization:
        
        The initialization of a Solution can happen on of these ways
            1. args is empty and an empty Solution is created
            2. args is a single State or list of States
            2. args is a single Grid or list of Grids
            3. args is a single Dimension or list of Dimensions
            4. args is a variable number of arguments that describes the 
               location of a file to be read in to initialize the object
        
        Input:
            - if args == () -> Empty Solution object
            - if args == States -> States are appended to states list
            - if args == Grids -> States are initialized with these Grids 
              and appended to states list
            - if args == Dimensions -> A single Grid with the given
              Dimensions is created, a state is initalized with this Grid
              and appended to the states list
            - if args == frame, format='ascii',path='./',file_prefix='fort'
    
    :Examples:

        >>> import pyclaw
        >>> x = pyclaw.Dimension('x',0.,1.,100)
        >>> grid = pyclaw.Grid((x))
        >>> state = pyclaw.State(grid,3,2)
        >>> solution = pyclaw.Solution(state)
    """

    # ========== Attributes ==================================================
    
    # ========== Properties ==================================================
    @property
    def state(self):
        r"""(:class:`State`) - Base state is returned"""
        return self.states[0]
    @property
    def grid(self):
        r"""(:class:`Grid`) - Base state's grid is returned"""
        return self.states[0].grid

    @property
    def t(self):
        r"""(float) - :attr:`State.t` of base state"""
        return self._get_base_state_attribute('t')
    @t.setter
    def t(self,value):
        self.set_all_states('t',value)

    @property
    def num_eqn(self):
        r"""(int) - :attr:`State.num_eqn` of base state"""
        return self._get_base_state_attribute('num_eqn')
    @property
    def mp(self):
        r"""(int) - :attr:`State.mp` of base state"""
        return self._get_base_state_attribute('mp')
    @property
    def mF(self):
        r"""(int) - :attr:`State.mF` of base state"""
        return self._get_base_state_attribute('mF')
    @property
    def q(self):
        r"""(ndarray(...,:attr:`State.num_eqn`)) - :attr:`State.q` of base state"""
        return self._get_base_state_attribute('q')
    @property
    def p(self):
        r"""(ndarray(...,:attr:`State.mp`)) - :attr:`State.p` of base state"""
        return self._get_base_state_attribute('p')
    @property
    def F(self):
        r"""(ndarray(...,:attr:`State.mF`)) - :attr:`State.F` of base 
                  state"""
        return self._get_base_state_attribute('F')

    @property
    def aux(self):
        r"""(ndarray(...,:attr:`State.num_aux`)) - :attr:`State.aux` of base 
                  state"""
        return self._get_base_state_attribute('aux')
    @aux.setter
    def aux(self, value): 
        if len(self.states) == 1: 
            setattr(self.states[0],'aux',value)

    @property
    def capa(self):
        r"""(ndarray(...)) - :attr:`State.capa` of base state"""
        return self._get_base_state_attribute('capa')
    @capa.setter
    def capa(self, value):
        if len(self.states) == 1: 
            setattr(self.states[0],'capa',value)

    @property
    def problem_data(self):
        r"""(dict) - :attr:`State.problem_data` of base state"""
        return self._get_base_state_attribute('problem_data')
    @problem_data.setter
    def problem_data(self, value):
        if len(self.states) == 1: 
            setattr(self.states[0],'problem_data',value)

    @property
    def num_aux(self):
        r"""(int) - :attr:`State.num_aux` of base state"""
        return self._get_base_state_attribute('num_aux')
    @property
    def ndim(self):
        r"""(int) - :attr:`Grid.ndim` of base state.grid"""
        return self._get_base_grid_attribute('ndim')
    @property
    def dimensions(self):
        r"""(list) - :attr:`Grid.dimensions` of base state.grid"""
        return self._get_base_grid_attribute('dimensions')
    @property
    def n(self):
        r"""(list) - :attr:`Grid.n` of base state.grid"""
        return self._get_base_grid_attribute('n')
    @property
    def name(self):
        r"""(list) - :attr:`Grid.name` of base state.grid"""
        return self._get_base_grid_attribute('name')
    @property
    def lower(self):
        r"""(list) - :attr:`Grid.lower` of base state.grid"""
        return self._get_base_grid_attribute('lower')
    @property
    def upper(self):
        r"""(list) - :attr:`Grid.upper` of base state.grid"""
        return self._get_base_grid_attribute('upper')
    @property
    def d(self):
        r"""(list) - :attr:`Grid.d` of base state.grid"""
        return self._get_base_grid_attribute('d')
    @property
    def units(self):
        r"""(list) - :attr:`Grid.units` of base state.grid"""
        return self._get_base_grid_attribute('units')
    @property
    def center(self):
        r"""(list) - :attr:`Grid.center` of base state.grid"""
        return self._get_base_grid_attribute('center')
    @property
    def edge(self):
        r"""(list) - :attr:`Grid.edge` of base state.grid"""
        return self._get_base_grid_attribute('edge')
    @property
    def p_center(self):
        r"""(list) - :attr:`Grid.p_center` of base state.grid"""
        return self._get_base_grid_attribute('p_center')
    @property
    def p_edge(self):
        r"""(list) - :attr:`Grid.p_edge` of base state.grid"""
        return self._get_base_grid_attribute('p_edge')
    @property
    def c_center(self):
        r"""(list) - :attr:`Grid.c_center` of base state.grid"""
        return self._get_base_grid_attribute('c_center')
    @property
    def c_edge(self):
        r"""(list) - :attr:`Grid.c_edge` of base state.grid"""
        return self._get_base_grid_attribute('c_edge')
        

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
                import inspect
                frame = arg[0]
                self.read(frame,**kargs)
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
        for state in self.states:
            if getattr(state,attr) is None or overwrite:
                setattr(state,attr,value) 
    
                    
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
    
    
    def __deepcopy__(self,memo={}):
        import copy
        # Create basic container
        result = self.__class__()
        result.__init__()
        
        # Populate the states
        for state in self.states:
            result.states.append(copy.deepcopy(state))
        
        return result
    
    
    # ========== IO Functions ================================================
    def write(self,frame,path='./',file_format='ascii',file_prefix=None,
                write_aux=False,options={},write_p=False):
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
        if isinstance(file_format,str):
            format_list = [file_format]
        elif isinstance(file_format,list):
            format_list = file_format
        if 'petsc' in format_list:
            from petclaw import io
        # Loop over list of formats requested
        for form in format_list:
            write_func = eval('io.write_%s' % form)
            if file_prefix is None:
                write_func(self,frame,path,write_aux=write_aux,
                            options=options,write_p=write_p)
            else:
                write_func(self,frame,path,file_prefix=file_prefix,
                                write_aux=write_aux,options=options,
                           write_p=write_p)
            msg = "Wrote out solution in format %s for time t=%s" % (form,self.t)
            logging.getLogger('io').info(msg)
        
    def read(self,frame,path='./_output',file_format='ascii',file_prefix=None,
                read_aux=True,options={}, **kargs):
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
        
        if file_format=='petsc':
            from petclaw import io
        path = os.path.expandvars(os.path.expanduser(path))
        read_func = eval('io.read_%s' % file_format)
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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
