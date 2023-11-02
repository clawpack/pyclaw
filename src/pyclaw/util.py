#!/usr/bin/env python
# encoding: utf-8
r"""
Pyclaw utility methods
"""

import time
import os
import sys
import subprocess
import logging
import tempfile
import inspect
import warnings

import numpy as np
import contextlib


LOGGING_LEVELS = {0:logging.CRITICAL,
                  1:logging.ERROR,
                  2:logging.WARNING,
                  3:logging.INFO,
                  4:logging.DEBUG}

def add_parent_doc(parent):
    """add parent documentation for a class""" 
    
    return """
    **Parent Class Documentation:**
    """ + parent.__doc__

def run_serialized(fun):
    """ Decorates a function to only run serially, even if called in parallel.

    In a parallel communicator, the first process will run while the remaining processes
    block on a barrier.  In a serial run, the function will be called directly.

    This currently assumes the global communicator is PETSc.COMM_WORLD, but is easily
    generalized.
    """

    try:
        from petsc4py import PETSc
        is_parallel = True
    except ImportError:
        is_parallel = False

    if is_parallel:
        rank = PETSc.COMM_WORLD.getRank()
        if rank == 0:
            def serial_fun(*args, **kwargs):
                fun(*args, **kwargs)
                PETSc.COMM_WORLD.Barrier()
        else:
            def serial_fun(*args, **kwargs):
                PETSc.COMM_WORLD.Barrier()
    else:
        def serial_fun(*args, **kwargs):
            fun(*args, **kwargs)

    return serial_fun

@run_serialized
def inplace_build(working_dir, warn=True):
    """Build missing extension modules with an in-place build.  This is a convenience
    function for PyClaw applications that rely on custom extension modules.  In its default
    mode, this function emits warnings indicating its actions.

    This function is safe to execute in parallel, it will only run on the 0 zero process while
    the other processes block.
    """

    if warn:
        warnings.warn("missing extension modules")
        warnings.warn("running python setup.py build_ext -i in %s" % working_dir)

    subprocess.check_call('python setup.py build_ext -i', shell=True, cwd=working_dir)

    if warn:
        warnings.warn("successfully executed python setup.py build_ext -i in %s" % working_dir)

def run_app_from_main(application,setplot=None):
    r"""
    Runs an application from pyclaw/examples/, automatically parsing command line keyword
    arguments (key=value) as parameters to the application, with positional
    arguments being passed to PETSc (if it is enabled).

    Perhaps we should take the PETSc approach of having a database of PyClaw
    options that can be queried for options on specific objects within the
    PyClaw runtime instead of front-loading everything through the application
    main...
    """

    # Arguments to the PyClaw should be keyword based, positional arguments
    # will be passed to PETSc
    petsc_args, pyclaw_kwargs = _info_from_argv(sys.argv)

    if 'use_petsc' in pyclaw_kwargs and pyclaw_kwargs['use_petsc']:
        import petsc4py
        petsc_args = [arg.replace('--','-') for arg in sys.argv[1:] if '=' not in arg]
        petsc4py.init(petsc_args)
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if sys.version_info >= (2, 7):
        app_kwargs = {key: value for key, value in pyclaw_kwargs.items() 
                      if not key in ('htmlplot','iplot')}
    else:
        # the above fails with Python < 2.7, so write it out...
        app_kwargs = {}
        for key,value in pyclaw_kwargs.items():
            if key not in ('htmlplot','iplot'):
                app_kwargs[key] = value

    claw=application(**app_kwargs)

    # Solve
    status = claw.run()

    # Plot results
    htmlplot = pyclaw_kwargs.get('htmlplot',False)
    iplot    = pyclaw_kwargs.get('iplot',False)
    outdir   = pyclaw_kwargs.get('outdir','./_output')
    if htmlplot:  
        if setplot is not None:
            pyclaw.plot.html_plot(outdir=outdir,setplot=setplot)
        else:
            pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     
        if setplot is not None:
            pyclaw.plot.interactive_plot(outdir=outdir,setplot=setplot)
        else:
            pyclaw.plot.interactive_plot(outdir=outdir)

    return claw

class VerifyError(Exception):
    pass

def gen_variants(application, verifier, kernel_languages=('Fortran',), disable_petsc=False, **kwargs):
    r"""
    Generator of runnable variants of a test application given a verifier

    Given an application, a script for verifying its output, and a
    list of kernel languages to try, generates all possible variants of the
    application to try by taking a product of the available kernel_languages and
    (petclaw/pyclaw).  For many applications, this will generate 4 variants:
    the product of the two main kernel languages ('Fortran' and 'Python'), against
    the the two parallel modes (petclaw and pyclaw).

    For more information on how the verifier function should be implemented,
    see util.test_app for a description, and util.check_diff for an example.

    All unrecognized keyword arguments are passed through to the application.
    """

    arg_dicts = build_variant_arg_dicts(kernel_languages, disable_petsc)

    for test_kwargs in arg_dicts:
        test_kwargs.update(kwargs)
        try:
            test_name = application.__module__
        except:
            test_name = inspect.getmodule(application)
        if 'solver_type' in test_kwargs:
            solver_info = 'solver_type={solver_type!s}, '
        else:
            solver_info = ''
        test = lambda: test_app(application, verifier, test_kwargs)
        test.description = '%s(%s)' % (test_name, str(kwargs))
        yield test
    return

def build_variant_arg_dicts(kernel_languages=('Fortran',), disable_petsc=False):
    import itertools

    # test petsc4py only if it is available
    use_petsc_opts = (False,)
    if not disable_petsc:
        try:
            import petsc4py
            use_petsc_opts=(True,False)
        except ImportError:
            pass # petsc starts disabled

    opt_names = 'use_petsc','kernel_language'
    opt_product = itertools.product(use_petsc_opts,kernel_languages)
    arg_dicts = [dict(list(zip(opt_names,argset))) for argset in opt_product]

    return arg_dicts

def test_app_variants(application, verifier, kernel_languages, **kwargs):

    arg_dicts = build_variant_arg_dicts(kernel_languages)

    for test_kwargs in arg_dicts:
        test_kwargs.update(kwargs)
        test_app(application, verifier, test_kwargs)
    return

def test_app(application, verifier, kwargs):
    r"""
    Test the output of a given application against its verifier method.

    This function performs the following two function calls::

        output = application(**kwargs)
        check_values = verifier(output)

    The verifier method should return None if the output is correct, otherwise
    it should return an indexed sequence of three items::

      0 - expected value
      1 - test value
      2 - string describing the tolerance type (abs/rel) and value.

    This information is used to present descriptive help if an error is detected.
    For an example verifier method, see util.check_diff

    """
    print(kwargs)

    if 'use_petsc' in kwargs and not kwargs['use_petsc']:
        try:
            # don't duplicate serial test runs
            from petsc4py import PETSc
            rank = PETSc.COMM_WORLD.getRank()
            if rank != 0:
                return
        except ImportError as e:
            pass

    claw = application(**kwargs)
    claw.run()
    check_values = verifier(claw)

    if check_values is not None:
        err = \
        """%s
********************************************************************************
verification function
%s
args                 : %s
norm of expected data: %s
norm of test data    : %s
test error           : %s
%s
********************************************************************************
""" % \
        (inspect.getsourcefile(application),
         inspect.getsource(verifier),
         kwargs,
         check_values[0],
         check_values[1],
         check_values[2],
         check_values[3])
        raise VerifyError(err)
    return

def check_diff(expected, test, **kwargs):
    r"""
    Checks the difference between expected and test values, return None if ok

    This function expects either the keyword argument 'abstol' or 'reltol'.
    """
    if 'delta' in kwargs:
        d = np.prod(kwargs['delta'])
    else:
        d = 1.0
    err_norm = d*np.linalg.norm(expected - test)
    expected_norm = d*np.linalg.norm(expected)
    test_norm = d*np.linalg.norm(test)
    if 'abstol' in kwargs:
        if err_norm < kwargs['abstol']: return None
        else: return (expected_norm, test_norm, err_norm,
                      'abstol  : %s' % kwargs['abstol'])
    elif 'reltol' in kwargs:
        if err_norm/expected_norm < kwargs['reltol']: return None
        else: return (expected_norm, test_norm, err_norm,
                      'reltol  : %s' % kwargs['reltol'])
    else:
        raise Exception('Incorrect use of check_diff verifier, specify tol!')

def check_solutions_are_same(sol_a,sol_b):
    assert len(sol_a.states) == len(sol_b.states)
    assert sol_a.t == sol_b.t
    for state in sol_a.states:
        for ref_state in sol_b.states:
            if ref_state.patch.patch_index == state.patch.patch_index:
                break

        # Required state attributes
        assert np.linalg.norm(state.q - ref_state.q) < 1.e-6 # Not sure why this can be so large
        if ref_state.aux is not None:
            assert np.linalg.norm(state.aux - ref_state.aux) < 1.e-16
        for attr in ['t', 'num_eqn', 'num_aux']:
            assert getattr(state,attr) == getattr(ref_state,attr)
        # Optional state attributes
        for attr in ['patch_index', 'level']:
            if hasattr(ref_state,attr):
                assert getattr(state,attr) == getattr(ref_state,attr)

        patch = state.patch
        ref_patch = ref_state.patch
        # Required patch attributes
        for attr in ['patch_index', 'level']:
            assert getattr(patch,attr) == getattr(ref_patch,attr)

        dims = patch.dimensions
        ref_dims = ref_patch.dimensions
        for dim, ref_dim in zip(dims,ref_dims):
            # Required dim attributes
            for attr in ['num_cells','lower','delta']:
                assert getattr(dim,attr) == getattr(ref_dim,attr)
            # Optional dim attributes
            for attr in ['units','on_lower_boundary','on_upper_boundary']:
                if hasattr(ref_dim,attr):
                    assert getattr(dim,attr) == getattr(ref_dim,attr)
# ============================================================================
#  F2PY Utility Functions
# ============================================================================
def compile_library(source_list,module_name,interface_functions=[],
                        local_path='./',library_path='./',f2py_flags='',
                        FC=None,FFLAGS=None,recompile=False,clean=False):
    r"""
    Compiles and wraps fortran source into a callable module in python.

    This function uses f2py to create an interface from python to the fortran
    sources in source_list.  The source_list can either be a list of names
    of source files in which case compile_library will search for the file in
    local_path and then in library_path.  If a path is given, the file will be
    checked to see if it exists, if not it will look for the file in the above
    resolution order.  If any source file is not found, an IOException is
    raised.

    The list interface_functions allows the user to specify which fortran
    functions are actually available to python.  The interface functions are
    assumed to be in the file with their name, i.e. claw1 is located in
    'claw1.f95' or 'claw1.f'.

    The interface from fortran may be different than the original function
    call in fortran so the user should make sure to check the automatically
    created doc string for the fortran module for proper use.

    Source files will not be recompiled if they have not been changed.

    One set of options of note is for enabling OpenMP, it requires the usual
    fortran flags but the OpenMP library also must be compiled in, this is
    done with the flag -lgomp.  The call to compile_library would then be:

    compile_library(src,module_name,f2py_flags='-lgomp',FFLAGS='-fopenmp')

    For complete optimization use:

    FFLAGS='-O3 -fopenmp -funroll-loops -finline-functions -fdefault-real-8'

    :Input:
     - *source_list* - (list of strings) List of source files, if these are
       just names of the source files, i.e. 'bc1.f' then they will be searched
       for in the default source resolution order, if an explicit path is
       given, i.e. './bc1.f', then the function will use that source if it can
       find it.
     - *module_name* - (string) Name of the resulting module
     - *interface_functions* - (list of strings) List of function names to
       provide access to, if empty, all functions are accessible to python.
       Defaults to [].
     - *local_path* - (string) The base path for source resolution, defaults
       to './'.
     - *library_path* - (string) The library path for source resolution,
       defaults to './'.
     - *f2py_flags* - (string) f2py flags to be passed
     - *FC* - (string) Override the environment variable FC and use it to
       compile, note that this does not replace the compiler that f2py uses,
       only the object file compilation (functions that do not have
       interfaces)
     - *FFLAGS* - (string) Override the environment variable FFLAGS and pass
       them to the fortran compiler
     - *recompile* - (bool) Force recompilation of the library, defaults to
       False
     - *clean* - (bool) Force a clean build of all source files
    """

    # Setup logger
    logger = logging.getLogger('f2py')
    temp_file = tempfile.TemporaryFile()
    logger.info('Compiling %s' % module_name)

    # Force recompile if the clean flag is set
    if clean:
        recompile = True

    # Expand local_path and library_path
    local_path = os.path.expandvars(local_path)
    local_path = os.path.expanduser(local_path)
    library_path = os.path.expandvars(library_path)
    library_path = os.path.expanduser(library_path)

    # Fetch environment variables we need for compilation
    if FC is None:
        if 'FC' in os.environ:
            FC = os.environ['FC']
        else:
            FC = 'gfortran'

    if FFLAGS is None:
        if 'FFLAGS' in os.environ:
            FFLAGS = os.environ['FFLAGS']
        else:
            FFLAGS = ''

    # Create the list of paths to sources
    path_list = []
    for source in source_list:
        # Check to see if the source looks like a path, i.e. it contains the
        # os.path.sep character
        if source.find(os.path.sep) >= 0:
            source = os.path.expandvars(source)
            source = os.path.expanduser(source)
            # This is a path, check to see if it's valid
            if os.path.exists(source):
                path_list.append(source)
                continue
            # Otherwise, take the last part of the path and try searching for
            # it in the resolution order
            source = os.path.split(source)

        # Search for the source file in local_path and then library_path
        if os.path.exists(os.path.join(local_path,source)):
            path_list.append(os.path.join(local_path,source))
            continue
        elif os.path.exists(os.path.join(library_path,source)):
            path_list.append(os.path.join(library_path,source))
            continue
        else:
            raise IOError('Could not find source file %s' % source)

    # Compile each of the source files if the object files are not present or
    # if the modification date of the source file is newer than the object
    # file's creation date
    object_list = []
    src_list = []
    for path in path_list:
        object_path = os.path.join(os.path.split(path)[0],
            '.'.join((os.path.split(path)[1].split('.')[:-1][0],'o')))

        # Check to see if this path contains one of the interface functions
        if os.path.split(path)[1].split('.')[:-1][0] in interface_functions:
            src_list.append(path)
            continue
        # If there are no interface functions specified, then all source files
        # must be included in the f2py call
        elif len(interface_functions) == 0:
            src_list.append(path)
            continue

        if os.path.exists(object_path) and not clean:
            # Check to see if the modification date of the source file is
            # greater than the object file
            if os.path.getmtime(object_path) > os.path.getmtime(path):
                object_list.append(object_path)
                continue
        # Compile the source file into the object file
        command = '%s %s -c %s -o %s' % (FC,FFLAGS,path,object_path)
        logger.debug(command)
        subprocess.call(command,shell=True,stdout=temp_file)
        object_list.append(object_path)

    # Check to see if recompile is needed
    if not recompile:
        module_path = os.path.join('.','.'.join((module_name,'so')))
        if os.path.exists(module_path):
            for src in src_list:
                if os.path.getmtime(module_path) < os.path.getmtime(src):
                    recompile = True
                    break
            for obj in object_list:
                if os.path.getmtime(module_path) < os.path.getmtime(obj):
                    recompile = True
                    break
        else:
            recompile = True

    if recompile:
        # Wrap the object files into a python module
        f2py_command = "f2py -c"
        # Add standard compiler flags
        f2py_command = ' '.join((f2py_command,f2py_flags))
        f2py_command = ' '.join((f2py_command,"--f90flags='%s'" % FFLAGS))
        # Add module names
        f2py_command = ' '.join((f2py_command,'-m %s' % module_name))
        # Add source files
        f2py_command = ' '.join((f2py_command,' '.join(src_list)))
        # Add object files
        f2py_command = ' '.join((f2py_command,' '.join(object_list)))
        # Add interface functions
        if len(interface_functions) > 0:
            f2py_command = ' '.join( (f2py_command,'only:') )
            for interface in interface_functions:
                f2py_command = ' '.join( (f2py_command,interface) )
            f2py_command = ''.join( (f2py_command,' :') )
        logger.debug(f2py_command)
        status = subprocess.call(f2py_command,shell=True,stdout=temp_file)
        if status == 0:
            logger.info("Module %s compiled" % module_name)
        else:
            logger.info("Module %s failed to compile with code %s" % (module_name,status))
            sys.exit(13)
    else:
        logger.info("Module %s is up to date." % module_name)

    temp_file.seek(0)
    logger.debug(temp_file.read())
    temp_file.close()

def construct_function_handle(path,function_name=None):
    r"""
    Constructs a function handle from the file at path.

    This function will attempt to construct a function handle from the python
    file at path.

    :Input:
     - *path* - (string) Path to the file containing the function
     - *function_name* - (string) Name of the function defined in the file
       that the handle will point to.  Defaults to the same name as the file
       without the extension.

    :Output:
     - (func) Function handle to the constructed function, None if this has
       failed.
    """
    # Determine the resulting function_name
    if function_name is None:
        function_name = path.split('/')[-1].split('.')[0]

    full_path = os.path.abspath(path)
    if os.path.exists(full_path):
        suffix = path.split('.')[-1]
        # This is a python file and we just need to read it and map it
        if suffix in ['py']:
            exec(compile(open(full_path).read(), full_path, 'exec'),globals())
            return eval('%s' % function_name)
        else:
            raise Exception("Invalid file type for function handle.")
    else:
        raise Exception("Invalid file path %s" % path)


#---------------------------------------------------------
def read_data_line(inputfile,num_entries=1,data_type=float):
#---------------------------------------------------------
    r"""
    Read data a single line from an input file

    Reads one line from an input file and returns an array of values

    inputfile: a file pointer to an open file object
    num_entries: number of entries that should be read, defaults to only 1
    type: Type of the values to be read in, they all must be the same type

    This function will return either a single value or an array of values
    depending on if num_entries > 1

    """
    l = []
    while  l==[]:  # skip over blank lines
        line = inputfile.readline()
        if line == '':
            raise IOError('*** Reached EOF in file %s' % inputfile.name)
        l = line.split()
    if num_entries == 1:  # This is a convenience for calling functions
        return data_type(l[0])
    val = np.empty(num_entries,data_type)
    if num_entries > len(l):
        print('Error in read_data_line: num_entries = ', num_entries)
        print('  is larger than length of l = ',l)
    val = [data_type(entry) for entry in l[:num_entries]]
    return val

#----------------------------------------
def convert_fort_double_to_float(number):
#----------------------------------------
    r"""
    Converts a fortran format double to a float

    Converts a fortran format double to a python float.

    number: is a string representation of the double.  Number should
    be of the form "1.0d0"

    """
    a = number.split('d')
    return float(a[0])*10**float(a[1])

#-----------------------------
def current_time(addtz=False):
#-----------------------------
    # determine current time and reformat:
    time1 = time.asctime()
    year = time1[-5:]
    day = time1[:-14]
    hour = time1[-13:-5]
    current_time = day + year + ' at ' + hour
    if addtz:
        current_time = current_time + ' ' + time.tzname[time.daylight]
    return current_time


def _method_info_from_argv(argv=None):
    """Command-line -> method call arg processing.

    - positional args:
            a b -> method('a', 'b')
    - intifying args:
            a 123 -> method('a', 123)
    - json loading args:
            a '["pi", 3.14, null]' -> method('a', ['pi', 3.14, None])
    - keyword args:
            a foo=bar -> method('a', foo='bar')
    - using more of the above
            1234 'extras=["r2"]'  -> method(1234, extras=["r2"])

    @param argv {list} Command line arg list. Defaults to `sys.argv`.
    @returns (<method-name>, <args>, <kwargs>)
    """
    import json
    if argv is None:
        argv = sys.argv

    method_name, arg_strs = argv[1], argv[2:]
    args = []
    kwargs = {}
    for s in arg_strs:
        if s.count('=') == 1:
            key, value = s.split('=', 1)
        else:
            key, value = None, s
        try:
            value = json.loads(value)
        except ValueError:
            pass
        if value=='True': value=True
        if value.lower()=='false': value=False
        if key:
            kwargs[key] = value
        else:
            args.append(value)
    return method_name, args, kwargs

def _info_from_argv(argv=None):
    """Command-line -> method call arg processing.

    - positional args:
            a b -> method('a', 'b')
    - intifying args:
            a 123 -> method('a', 123)
    - json loading args:
            a '["pi", 3.14, null]' -> method('a', ['pi', 3.14, None])
    - keyword args:
            a foo=bar -> method('a', foo='bar')
    - using more of the above
            1234 'extras=["r2"]'  -> method(1234, extras=["r2"])

    @param argv {list} Command line arg list. Defaults to `sys.argv`.
    @returns (<method-name>, <args>, <kwargs>)
    """
    import json
    if argv is None:
        argv = sys.argv

    arg_strs = argv[1:]
    args = []
    kwargs = {}
    for s in arg_strs:
        if s.count('=') == 1:
            key, value = s.split('=', 1)
        else:
            key, value = None, s
        try:
            value = json.loads(value)
        except ValueError:
            pass
        if value=='True': value=True
        if value=='False': value=False
        if key:
            kwargs[key] = value
        else:
            args.append(value)
    return args, kwargs

def _arguments_str_from_dictionary(options):
    """
    Convert method options passed as a dictionary to a str object
    having the form of the method arguments
    """
    option_string = ""
    for k in options:
        if isinstance(options[k], str):
            option_string += k+"='"+str(options[k])+"',"
        else:
            option_string += k+"="+str(options[k])+","
    option_string = option_string.strip(',')

    return option_string
#-----------------------------
class FrameCounter:
#-----------------------------
    r"""
    Simple frame counter

    Simple frame counter to keep track of current frame number.  This can
    also be used to keep multiple runs frames separated by having multiple
    counters at once.

    Initializes to 0
    """
    def __init__(self):
        self.__frame = 0

    def __repr__(self):
        return str(self.__frame)

    def increment(self):
        r"""
        Increment the counter by one
        """
        self.__frame += 1
    def set_counter(self,new_frame_num):
        r"""
        Set the counter to new_frame_num
        """
        self.__frame = new_frame_num
    def get_counter(self):
        r"""
        Get the current frame number
        """
        return self.__frame
    def reset_counter(self):
        r"""
        Reset the counter to 0
        """
        self.__frame = 0
