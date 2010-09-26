

.. _sharing:

##########################
Saving and sharing results
##########################

Clawpack now includes some tools to help facilitate archiving and sharing
results that you have obtained with this software.  
These make it relatively easy to generate
a set of webpages such as those seen when browsing the examples listed at
:ref:`apps` or :ref:`book`.

These webpages can easily be posted on your own website to be viewed by
others if you wish, for example to share on-going work with collaborators or
to supplement a journal article.

Making html files from codes and README files
=============================================

In a directory with a standard Clawpack Makefile the command::

   $ make .htmls

will create html files from all the .f, .py, .data files in this directory
and will also convert README.txt into README.html.   This uses the python
script `$CLAW/doc/clawcode2html.py <claw/doc/clawcode2html.py>`_
to convert each file.  

This script is based on a more general script 
`mathcode2html <http://www.amath.washington.edu/~rjl/mathcode2html/>`_ and
the documentation of that contains some examples of formating.  



Making webpages of plots
========================

The "make .plots" option available via the standard Makefiles will create a
set of webpages illustrating the plots and allowing easy navigation between
frames.  These webpages also allow viewing all frames of a plot as an
animation (via javascript within the browser).  For example see

   * `$CLAW/apps/advection/2d/example1/_plots/_PlotIndex.html
     <claw/apps/advection/2d/example1/_plots/_PlotIndex.html>`_

See :ref:`plotting` for more details on how to specify what plots will
appear on these webpages.

Sharing your results
====================

To make it easy for others to view your code and the resulting plots, you
can simply copy the example directory (containg the code, the .html files,
and the _plots subdirectory)  to your publicly visible web pages.  

