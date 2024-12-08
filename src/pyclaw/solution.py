#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing all Pyclaw solution objects
"""

import os
import logging

from .geometry import Patch, Dimension, Domain

# ============================================================================
#  Solution Class
# ============================================================================
class Solution(object):
    r"""
    Pyclaw patch container class
        
    :Input and Output:
    
        Input and output of solution objects is handle via the io package.
        Solution contains the generic methods :meth:`write`, :meth:`read` and
        :meth:`plot` which then figure out the correct method to call.  Please
        see the io package for the particulars of each format and method and 
        the methods in this class for general input and output information.
    
    :Properties:
    
        If there is only one state and patch belonging to this solution, 
        the solution will appear to have many of the attributes assigned to its
        one state and patch.  Some parameters that have in the past been
        parameters for all patch,s are also reachable although Solution does not
        check to see if these parameters are truly universal.

        Patch Attributes:
            'dimensions'
        State Attributes:
            't','num_eqn','q','aux','capa','problem_data'
            
            
    :Initialization:
        
        The initialization of a Solution can happen one of these ways
        
            1. `args` is empty and an empty Solution is created
            2. `args` is an integer (the number of components of q), a single
               State, or a list of States and is followed
               by the appropriate :ref:`geometry <pyclaw_geometry>` object
               which can be one of:
                
                 - (:class:`~pyclaw.geometry.Domain`)
                 - (:class:`~pyclaw.geometry.Patch`) - A domain is created
                   with the patch or list of patches provided.
                 - (:class:`~pyclaw.geometry.Dimension`) - A domain and 
                   patch is created with the dimensions or list of 
                   dimensions provided.
            3. `args` is a variable number of keyword arguments that describes the
               location of a file to be read in to initialize the object
    
    :Examples:

        >>> import clawpack.pyclaw as pyclaw
        >>> x = pyclaw.Dimension('x',0.,1.,100)
        >>> domain = pyclaw.Domain((x))
        >>> state = pyclaw.State(domain,3,2)
        >>> solution = pyclaw.Solution(state,domain)
    """
    def __getattr__(self, key):
        if key in ('t','num_eqn','mp','mF','q','p','F','aux','capa',
                   'problem_data','num_aux',
                   'num_dim', 'p_centers', 'p_edges', 'c_centers', 'c_edges',
                   'num_cells', 'lower', 'upper', 'delta', 'centers', 'edges',
                   'gauges', 'num_eqn', 'num_aux', 'grid', 'problem_data'):
            return self._get_base_state_attribute(key)
        else:
            raise AttributeError("'Solution' object has no attribute '"+key+"'")

    def __setattr__(self, key, value):
        if key in ('t','mp','mF'):
            self.set_all_states(key,value)
        else:
            self.__dict__[key] = value

    # ========== Attributes ==================================================
    
    # ========== Properties ==================================================
    @property
    def state(self):
        r"""(:class:`State`) - Base state is returned"""
        return self.states[0]
    @property
    def patch(self):
        r"""(:class:`Patch`) - Base state's patch is returned"""
        return self.domain.patch

    @property
    def start_frame(self):
        r"""(int) - : Solution start frame number in case the `Solution`
        object is initialized by loading frame from file"""
        return self._start_frame
    _start_frame = 0


    # ========== Class Methods ===============================================
    def __init__(self,*arg,**kargs):
        r"""Solution Initiatlization Routine
        
        See :class:`Solution` for more info.
        """

        # select package to build solver objects from, by default this will be
        # the package that contains the module implementing the derived class
        # for example, if Solution is implemented in 'clawpack.petclaw.solution', then 
        # the computed claw_package will be 'clawpack.petclaw'

        import sys
        if 'claw_package' in list(kargs.keys()):
            claw_package = kargs['claw_package']
        else:
            claw_package = None

        if claw_package is not None and claw_package in sys.modules:
            self.claw_package = sys.modules[claw_package]
        else:
            def get_clawpack_dot_xxx(modname): return modname.rpartition('.')[0]
            claw_package_name = get_clawpack_dot_xxx(self.__module__)
            if claw_package_name in sys.modules:
                self.claw_package = sys.modules[claw_package_name]
            else:
                raise NotImplementedError("Unable to determine solver package, please provide one")

        State = self.claw_package.State

        self.states = []
        self.domain = None
        if len(arg) == 1:
            # Load frame
            frame = int(arg[0])
            if ('count_from_zero' in kargs):
                if (kargs['count_from_zero'] == True):
                    self._start_frame = 0
                else:
                    self._start_frame = frame
                kargs.pop('count_from_zero')
            self.read(frame,**kargs)

        elif len(arg) == 2:
            #Set domain
            if isinstance(arg[1],Domain):
                self.domain = arg[1]
            else:
                if not isinstance(arg[1],(list,tuple)):
                    arg[1] = list(arg[1])
                if isinstance(arg[1][0],Dimension):
                    self.domain = Domain(Patch(arg[1]))
                elif isinstance(arg[1][0],Patch):
                    self.domain = Domain(arg[1])
                else:
                    raise Exception("Invalid arguments for Solution initialization.")

            #Set state
            if isinstance(arg[0],State):
                # Single State
                self.states.append(arg[0])
            elif isinstance(arg[0],(list,tuple)):
                if isinstance(arg[0][0],State):
                    # List of States
                    self.states = arg[0]
                elif isinstance(arg[0][0],int):
                    self.states = State(self.domain,arg[0][0],arg[0][1])
                else:
                    raise Exception("Invalid arguments for Solution initialization.")
            elif isinstance(arg[0],int):
                self.states.append(State(self.domain,arg[0]))
            if self.states == [] or self.domain is None:
                raise Exception("Invalid arguments for Solution initialization.")
        elif len(arg) == 0:
            if 'frame' in kargs:
                frame = int(kargs.pop('frame'))
                self.read(frame,**kargs)
            elif not kargs:
                pass  # With no arguments, initialize empty solution
            else:
                raise Exception("Invalid arguments for Solution initialization.")
        else:
            raise Exception("Invalid arguments for Solution initialization.")
                
                
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
        result.domain = copy.deepcopy(self.domain)
        
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
         - *write_aux* - (book) Write the auxiliary array out as well if 
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
                print("directory already exists, ignoring")  

        # Call the correct write function based on the output format
        if isinstance(file_format,str):
            format_list = [file_format]
        elif isinstance(file_format,list):
            format_list = file_format


        # Loop over list of formats requested
        for form in format_list:
            write_func = self.get_write_func(form)

            if file_prefix is None:
                write_func(self,frame,path,write_aux=write_aux,
                            options=options,write_p=write_p)
            else:
                write_func(self,frame,path,file_prefix=file_prefix,
                                write_aux=write_aux,options=options,
                           write_p=write_p)
            msg = "Wrote out solution in format %s for time t=%s" % (form,self.t)
            logging.getLogger('pyclaw.fileio').info(msg)

        
    def read(self, frame, path='./_output', file_format=None, 
                          file_prefix='fort', read_aux=True, options={}, **kargs):
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
           ``default = './_output'``
         - *file_format* - (string) Format of the file, should match on of the 
           modules inside of the io package.  ``default = None``
           but now attempts to read from header file (as of v5.9.0).
         - *file_prefix* - (string) Name prefix in front of all the files, 
           defaults to whatever the format defaults to, e.g. fort for ascii
         - *options* - (dict) Dictionary of optional arguments dependent on 
           the format being read in.  ``default = {}``
            
        :Output:
         - (bool) - True if read was successful, False otherwise
        """

        from clawpack.pyclaw.fileio.ascii import read_t
        file_format2 = None
        
        try:
            [t,num_eqn,nstates,num_aux,num_dim,num_ghost,file_format2] = \
                 read_t(frame,path,file_prefix=file_prefix)
        except:
            pass

        if file_format2 is not None:
            # value was read in from file, use it:
            file_format = file_format2

        read_func = self.get_read_func(file_format)

        if file_format == 'petsc':
            options['format'] = 'binary'
        else:
            options['format'] = file_format

        path = os.path.expandvars(os.path.expanduser(path))
        if file_prefix is None:
            read_func(self,frame,path,read_aux=read_aux,options=options)
        else:
            read_func(self,frame,path,file_prefix=file_prefix,
                                    read_aux=read_aux,options=options)
        logging.getLogger('pyclaw.fileio').info("Read in solution for time t=%s" % self.t)


    def get_read_func(self, file_format):
        if file_format[:6] == 'binary':
            # could be 'binary64' or 'binary32'
            import clawpack.pyclaw.fileio.binary
            return clawpack.pyclaw.fileio.binary.read
        elif file_format == 'ascii':
            import clawpack.pyclaw.fileio.ascii
            return clawpack.pyclaw.fileio.ascii.read
        elif file_format in ('hdf','hdf5'):
            import clawpack.pyclaw.fileio.hdf5
            return clawpack.pyclaw.fileio.hdf5.read
        elif file_format == 'petsc':
            try:
                import clawpack.petclaw.fileio
                return clawpack.petclaw.fileio.petsc.read
            except AttributeError as e:
                try:
                    from petsc4py import PETSc
                except ImportError:
                    raise ImportError("petsc4py is required for reading petsc format files, but petsc4py could not be imported.")
                raise e
        elif file_format == 'forestclaw':
            import clawpack.forestclaw.fileio.ascii
            return clawpack.forestclaw.fileio.ascii.read
        else:
            raise ValueError("File format %s not supported." % file_format)


    def get_write_func(self, file_format):
        if file_format == "forestclaw":
            import clawpack.forestclaw.fileio.ascii
            return clawpack.forestclaw.fileio.ascii.write
        elif file_format == "vtk":
            import clawpack.pyclaw.fileio.claw_vtk
            return clawpack.pyclaw.fileio.claw_vtk.write
        else:
            try:
                import importlib
                return importlib.import_module("clawpack.pyclaw.fileio.%s"
                                               % file_format).write
            except:
                raise ValueError("File format %s not found." % file_format)

        
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
