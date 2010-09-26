
.. _makefiles:


*************************************
Clawpack Makefiles
*************************************

In most directories with a `Makefile` you can type::

    $ make help

to find out what options are available.

Applications directory Makefiles
--------------------------------

In applications directories, compiling and running the code can usually be
accomplished via::

    $ make .output

This checks dependencies using the data of the hidden file `.output` that is
created after the code has successfully run.  If any Fortran codes have been
modified since this date, the code is first recompiled.  If the `setrun.py`
script has been changed more recently, then the data files are first
recreated.

Variables
+++++++++

A number of variables are defined in the Makefiles of application
directories.  For example, output is directed to the subdirectory specified
by the variable `OUTDIR`.  To change this, simply modify the Makefile before
typing "make .output".  Alternatively, you can modify the variable from the
command line, e.g.::

    $ make .output OUTDIR=run1

to direct output to a subdirectory named `run1`.

Compiler flags
++++++++++++++

Compiler flags can be changed by modifying the `FFLAGS` variable in the
Makefile.  If you change compiler flags you will generally need to recompile
all the Fortran files and the Makefile dependencies will not detect this.
To force recompilation of all files, use the "make new" option, e.g. to
recompile with the `-g` flag for debugging::

    $ make new FFLAGS=-g



.. warning::
   Some significant changes to Makefiles are contemplated for Clawpack 5.0.



