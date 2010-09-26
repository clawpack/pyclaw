
.. _claw43to44:

##########################################
Converting from Clawpack 4.3 to 4.4 or 4.5
##########################################


The minimal change needed to make a 4.3 application run in 4.4 (or 4.5)
is to rename
the *clawNez.data* file as *claw.data* and add a line to the top with a
single value *N*, the number of space dimensions.

To take full advantage of new features available in Clawpack 4.4 you may
also want to:

 * Create a *setrun.py* file to automatically generate *claw.data*,

 * Modify *setrun.py* to also create  *setprob.data* or other data files you
   read in. In this case you must alse modify *setprob.f* to use the 
   library routine *opendatafile* instead of the f77 function *open* as 
   illustrated in examples.  This routine skips over the warning message
   generated at the top of data files.

 * Switch to Python graphics by creating an appropriate *setplot.py* file,

 * Rewrite the Makefile to give options 'make .output', 'make .plots', etc.

 * Create a README.txt file that has links to the various other files in
   this directory, that will be converted to README.html by 'make .htmls'

Python conversion tool
----------------------

A first pass at all of the conversions listed above can often be achieved by
typing::

    $ python $CLAW/python/convert43to44.py

in your application directory.  You should then inspect the files generated
and fix any broken links, etc.

Currently this works for 1d and 2d applications on a single grid, 
but not for *amrclaw*.
