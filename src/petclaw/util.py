#!/usr/bin/env python
# encoding: utf-8
r"""
petclaw utility methods

:Authors:
    Kyle T. Mandli (2008-08-07) Initial version
    
    Randall J. LeVeque (2009-01-01) Added svn revision
    
    Kyle T. Mandli (2009-03-01) Added topo file utilities
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#      Copyright (C) 2009 Randall J. LeVeque <rjl@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import time
import os,sys,shutil,glob
import re
import subprocess
import logging
import tempfile

import numpy as np

# ============================================================================
#  Geoclaw Topography Utility Functions
# ============================================================================
def write_topo_file(path,z_func,nx,ny,lower_coord,upper_coord,topo_type=1):
    r"""
    Write out a topo1 file to the file at path.
    
    This utility function is for use with GeoClaw.  Please to refer to 
    GeoClaw's documentation for more details.
    
    :Input:
     - *path* - (string) Path to file
     - *z_func* - (func) Funtion determining z values (depth or topography)
     - *nx* - (int) Number of x points
     - *ny* - (int) Number of y points
     - *lower_coord* - (float) Lower left coordinate of the domain
     - *upper_coord* - (float) Upper right coordinate of the domain
     - *topo_type* - (int) Type of topography file to write, can be one of the
       following types.
       
    +---------------+--------------------------------------------------------+
    | Topo Type     | Description                                            |
    +===============+========================================================+
    | 1             | Standard GIS format:  3 columns with longitude,        |
    |               | latitude, and z in meters or x, y, z in one row        |
    +---------------+--------------------------------------------------------+
    | 2             | One z value per line with header                       |
    +---------------+--------------------------------------------------------+
    | 3             | One row of z values per line with header               |
    +---------------+--------------------------------------------------------+

    Header format is ::
    
        xxxx     ncols
        xxxx     nrows
        xxxx     xllcorner
        xxxx     yllcorner
        xxxx     cellsize
        xxxx     NODATA_value


    z values are positive for topography, negative for bathymetry.

    Values are ordered from upper left corner, moving east across each row
    (y fixed, x increasing) and then south (y decreasing).

    The region covered by this file should be refined to at least minlevel
    and at most maxlevel over the time interval specified.

    If a point is covered by several regions, the largest values of minlevel 
    and maxlevel from all such regions are used.

    Refinement may also be forced or allowed by the minlevel and maxlevel
    parameters for the regions specified in setregions.data or 
    setmovetopo.data.

    At points that lie far from shore (where h is greater than depthdeep),
    refinement is only allowed to level maxleveldeep.  These parameters are 
    set in setgeo.data.

    In setting topography in a computation, the finest topography available
    will be used for each grid cell.  If a grid cells lies partially but
    not entirely in a topo file then a coarser topo file will be used for
    the remainder of the cell (with area-weighted averaging).

    The order of topo files in the list above should not matter.
    """
    
    # Open file
    file = open(path,'w')
    
    # Calculate header values
    nrows = ny
    ncols = nx
    dx = (upper_coord[0]-lower_coord[0]) / nx
    dy = (upper_coord[1]-lower_coord[1]) / ny
    
    if not(dx == dy):
        raise ValueError("dx != dy, dx = %s, dy = %s" % (dx,dy))
        return 13
    cellsize = dx
    
    # Write header
    header = '%s ncols' % ncols
    header = '\n'.join((header,'%s nrows' % nrows))
    header = '\n'.join((header,'%s xllcorner' % lower_coord[0]))
    header = '\n'.join((header,'%s yllcorner' % lower_coord[1]))
    header = '\n'.join((header,'%s cellsize' % cellsize))
    header = '\n'.join((header,'%s NODATA_value\n' % -9999.0))
    file.write(header)
    
    # Write out function give topography type requested
    if topo_type == 1:
        for y in -np.linspace(-upper_coord[1],-lower_coord[1],ny):
            for x in np.linspace(lower_coord[0],upper_coord[0],nx):
                file.write("%s %s %s\n" % (x,y,z_func(x,y)))
    elif topo_type == 2:
        for y in -np.linspace(-upper_coord[1],-lower_coord[1],ny):
            for x in np.linspace(lower_coord[0],upper_coord[0],nx):
                file.write("%s\n" % z_func(x,y))
    elif topo_type == 3:
        for y in -np.linspace(-upper_coord[1],-lower_coord[1],ny):
            for x in np.linspace(lower_coord[0],upper_coord[0],nx):
                file.write("%s " % z_func(x,y))
            file.write("\n")
     
    # Close file
    file.close()
    
def read_topo_file(path,topo_type,verbose=False):
    r"""
    Read in a topography file from path of type topo_type.
    
    This utility function is for use with GeoClaw.  Please to refer to 
    GeoClaw's documentation for more details.
    
    :Input:
     - *path* - (string) Path to topography file
     - *topo_type* - (int) Type of topography file, see 
       :meth:`write_topo_file` for valid types.
     - *verbose* - (bool) Verbose output, ``default = False``
    
    :Output:
     - (ndarray(:)) X-Coordinate array
     - (ndarray(:)) Y-Coordinate array
     - (ndarray(:,:)) Z values stored in topography file
    """
    
    topo_file = open(path,'r')
    
    # Read in topography file header
    nx = int(topo_file.readline().split()[0])
    ny = int(topo_file.readline().split()[0])
    xll = float(topo_file.readline().split()[0])
    yll = float(topo_file.readline().split()[0])
    dx = float(topo_file.readline().split()[0])
    nodata_value = float(topo_file.readline().split()[0])
    
    if verbose:
        print "Header:"
        print "  (nx,ny) = (%s,%s)" % (nx,ny)
        print "  (xl,yl) = (%s,%s)" % (xll,yll)
        print "  (dx,dy) = (%s,%s)" % (dx,dx)
        print "  nodata  = %s" % nodata_value
    
    # Create Z matrix
    dy = dx
    x = np.linspace(xll,xll + nx * dx,nx)
    y = np.linspace(yll,yll + ny * dy,ny)
    Z = np.empty((nx,ny))
    
    if topo_type == 1:
        for j in -np.linspace(-ny+1,0,ny):
            for i in xrange(0,nx):
                temp = topo_file.readline()
                Z[i,j] = float(temp.split()[-1])
    elif topo_type == 2:
        for j in -np.linspace(-ny+1,0,ny):
            for i in xrange(0,nx):
                Z[i,j] = float(topo_file.readline())
    elif topo_type == 3:
        for j in -np.linspace(-ny+1,0,ny):
            line = topo_file.readline().split()
            for i in xrange(0,nx):
                Z[i,j] = float(line[i])
    else:
        print >> sys.stderr, "Invalid topo type %s" % topo_type
        return None
    
    return x,y,Z


def create_topo_func(loc,verbose=False):
    """Given a set of (x,z) locations, create a lambda function
    
    Create a lambda function that when evaluated will give the topgraphy 
    height at the point (x,y).
    
    :Example: 
    >>> f = create_topo_profile_func(loc)
    >>> b = f(x,y)
    
    :Input:
     - *loc* (list) - Create a topography file with the profile denoted by the
       tuples inside of loc.  A sample set of points are shown below.  Note 
       that the first value of the list is the x location and the second is 
       the height of the topography.
        
        z (m)
        ^                                                  o loc[5]  o
        |                                                    
        |                                          loc[4]   
        |--------------------------------------------o-----> x (m) (sea level)
        |                                            
        |                                o loc[2] o loc[3]
        |                         
        |                         
        |                           o loc[1]
        |           
        |                               
        |__________________o loc[0]
        0.0               
        
        
    """
    
    cmd_str = "lambda x,y: (x < %s) * %s" % (loc[0][0],loc[0][1])
    for i in xrange(0,len(loc)-1):
        loc_str = " + (%s < x <= %s)" % (loc[i][0],loc[i+1][0])
        loc_str = "".join((loc_str," * ((%s - %s) " % (loc[i][1],loc[i+1][1])))
        loc_str = "".join((loc_str," / (%s - %s)" % (loc[i][0],loc[i+1][0])))
        loc_str = "".join((loc_str," * (x - %s) + %s)" % (loc[i][0],loc[i][1])))
        cmd_str = "".join((cmd_str,loc_str))
    cmd_str = "".join((cmd_str," + (%s < x) * %s" % (loc[-1][0],loc[-1][1])))
    
    if verbose:
        print cmd_str
    return eval(cmd_str)
    

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
        if os.environ.has_key('FC'):
            FC = os.environ['FC']
        else:
            FC = 'gfortran'
      
    if FFLAGS is None:          
        if os.environ.has_key('FFLAGS'):
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
            execfile(full_path,globals())
            return eval('%s' % function_name)
        else:
            raise Exception("Invalid file type for function handle.")
    else:
        raise Exception("Invalid file path %s" % path)

#-----------------------------
class FrameCounter:
#-----------------------------
    r"""
    Simple frame counter

    Simple frame counter to keep track of current frame number.  This can
    also be used to keep multiple runs frames seperated by having multiple 
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


#---------------------------------------------------------
def read_data_line(inputfile,num_entries=1,type='float'):
#---------------------------------------------------------
    r"""
    Read data a single line from an input file

    Reads one line from an input file and returns an array of values

    inputfile: a file pointer to an open file object
    num_entries: number of entries that should be read, defaults to only 1
    type: Type of the values to be read in, they all most be the same type

    This function will return either a single value or an array of values
    depending on if num_entries > 1

    """
    l = []
    while  l==[]:  # skip over blank lines
        line = inputfile.readline()
        l = line.split()
    val = np.empty(num_entries,type)
    if num_entries > len(l):
        print 'Error in read_data_line: num_entries = ', num_entries
        print '  is larger than length of l = ',l
    try:
        for i in range(num_entries):
            exec("val[i] = %s(l[i])" % type)
        if num_entries == 1:  # This is a convenience for calling functions    
            return val[0]
        return val
    except(ValueError):
        print "Invalid type for the %s value in %s" % (i,l)
        print "  type = ",type
        return None
    except:
        raise

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
    print float(a[0])*10**float(a[1])
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



#---------------------------------
def svn_revision(dir="CLAW"):
#---------------------------------
    r"""
    Determine the svn revision number of the version of code being used.
    If dir="CLAW" it returns the revision number of the claw directory.
    This only checks the top most level and assumes this is accurate.
    """
   
    if dir=="CLAW":
        dir  = os.environ['CLAW']
    svn_entries = dir + '/.svn/entries'
    try:
        f = open(svn_entries,'r')
        lines = f.readlines()
        revision = int(lines[3])  # fourth line of file, converted to int
        f.close()
    except:
        revision = None
    
    return revision
