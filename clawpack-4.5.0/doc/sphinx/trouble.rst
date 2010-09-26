
.. _trouble:


*************************************
Troubleshooting
*************************************

Installation
++++++++++++

.. _trouble_makeexe:

Trouble running "make .exe"
---------------------------

If the code does not compile, check the following:

 * Make sure your environment variable `CLAW` is set properly::

    $ printenv CLAW

   to print the value.  
   The Makefiles use this variable to find the common Makefile and
   library routines.

 * Make sure your environment variable `FC` is set properly.  This
   should be set to
   the command used to invoke the Fortran compiler.  In many Makefiles
   this is set to `gfortran` by default if the user has not set it,
   via a line of the form::

     FC ?= gfortran   

   but this is ignored if the variable has been set by the user.
   If you get an error like::

    make[1]: gfortran: No such file or directory

   then the gfortran compiler is not being found.



.. _trouble_makedata:

Trouble running "make .data"
------------------------------

If you get the Python error::

    ImportError: No module named pyclaw

it's possible that your environment variable `PYTHONPATH` does not
include $CLAW/python in the path.  If you followed the instructions
at :ref:`setenv` then this variable should be set properly.   Recall
that this has to be set in any new shell you use for Clawpack.

If there are errors in the `setrun` function (usually defined in
`setrun.py`) then the these may show up when you try to "make .data"
since this function must be executed.


.. _trouble_makeoutput:

Trouble running "make .output"
------------------------------

If you want to re-run the code and you get::

    $ make .output
    make: `.output' is up to date.

then you can force it to run again by removing the file `.output`::

    $ rm -f .output
    $ make .output

This happens for example if you changed something that you know
will affect the output but that isn't in the Makefile's set of
dependencies.



.. _trouble_makeplots:

Trouble running "make .plots"
------------------------------
   
The Python plotting routines require `NumPy` and `matplotlib`.  See 
:ref:`python` for information on installing this.

If there are errors in the `setplot` function (usually defined in
`setplot.py`) then the these may show up when you try to "make .output"
since this function must be executed.


