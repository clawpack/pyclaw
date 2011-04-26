  .. _pyclaw_io:

***************************
Pyclaw Input/Output Package
***************************

Pyclaw supports the following input and output formats:

 * ASCII_ - ASCII file I/O, supports traditional clawpack format files
 * HDF5_ - HDF5 file I/O
 * NetCDF_ - NetCDF file I/O, support for NetCDF3 and NetCDF4 files
 
Each module contains two main routines ``read_<format>`` and 
``write_<format>`` which :class:`~pyclaw.Solution` can call with the 
appropriate <format>.  In order to create a new file I/O extension the calling
signature must match

::

    read_<format>(solution,frame,path,file_prefix,write_aux,options)
    
where the the inputs are
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be 
       output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name.
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out.    
     - *options* - (dict) Optional argument dictionary

and

::

    write_<format>(solution,frame,path,file_prefix,write_aux,options)

where the inputs are
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be 
       output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name.
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out.   
     - *options* - (dict) Optional argument dictionary.
     
Note that both allow for an ``options`` dictionary that is format specific
and should be documented thoroughly.  For examples of this usage, see the 
HDF5_ and NetCDF_ modules.
 
HDF5_ and NetCDF_ support require installed libraries in order to work, please
see the respective modules for details on how to obtain and install the
libraries needed.

  .. note::
  
     Pyclaw automatically detects the availability of HDF5 and NetCDF file 
     support and will warn you if you try and use them without the proper
     libraries.
     
.. _ASCII:

:mod:`pyclaw.io.ascii`
======================

.. automodule:: pyclaw.io.ascii
    :members:
    
.. _HDF5:

:mod:`pyclaw.io.hdf5`
=====================

.. automodule:: pyclaw.io.hdf5
    :members:

.. _NetCDF:

:mod:`pyclaw.io.netcdf`
=======================

.. automodule:: pyclaw.io.netcdf
    :members: