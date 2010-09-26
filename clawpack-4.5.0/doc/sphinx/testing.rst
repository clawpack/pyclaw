
.. _testing:

===================================================================
Regression tests
===================================================================

We are working on a full suite of regression tests for Clawpack.

For now, the best way to test that everything works properly is to use the
scripts `make_all.py`, `make_gallery.py`, and `make_htmls.py` found in the
directory `$CLAW/python <claw/python/README.html>`_.   These will run all
the examples in the `$CLAW/apps <claw/apps>`_ and `$CLAW/book <claw/book>`_
directories, make plots of the results, and a create a
gallery of thumbnails that should look like those in the 
`online Applications Galleries
<http://kingkong.amath.washington.edu/clawpack/users/apps.html>`_

Creating the gallery
--------------------


 * Make sure your envirornment variables are set properly, see
   :ref:`setenv`.

 * Also set the environment variable `CLAW_TOPO_DOWNLOAD` to True to avoid
   being prompted before topography files are downloaded.

 * Start the python webserver to best view html files of results.
   See :ref:`startserver`.

 * First test your installation by running a few examples by hand, for
   example::

    $ cd $CLAW/apps/advection/1d
    $ make .plots

   See :ref:`first_test` for more details.

 * Next test that the make_all.py script works by running a 
   few examples, e.g. those for 1d advection and acoustics::

    $ python $CLAW/python/make_all.py   $CLAW/apps/advection/1d
    $ python $CLAW/python/make_all.py   $CLAW/apps/acoustics/1d

 * Compile some files in the libraries needed to make sure modules are available.
   This step should be eliminated but sometimes it's needed... ::

    $ cd $CLAW/amrclaw/2d/lib
    $ make
    $ cd $CLAW/geoclaw/2d/lib
    $ make

 * To run all examples and make plots::
    
     $ cd $CLAW
     $ python $CLAW/python/make_all.py   

   .. warning:: Running all the examples may take an hour or two depending on your
      computer.
    
 * Also make all the .html files for README's and html versions of source
   files::

     $ python $CLAW/python/make_htmls.py   

 * Create galleries::

     $ python $CLAW/python/make_gallery.py

The galleries should then be in $CLAW/myclaw/gallery.
If you have your web server running (see :ref:`startserver`), they will be at
`<http://localhost:50005/doc/gallery>`_.
