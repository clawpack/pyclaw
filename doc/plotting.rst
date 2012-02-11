==========
Plotting
==========
PyClaw currently supports two output formats: ASCII and a custom binary format
referred to internally as `format=petsc`.  PyClaw relies on the 
`VisClaw package <http://github.com/clawpack/visclaw/>`_ for easy plotting, although
it is of course possible to load the output into other visualization packages.
VisClaw supports 1D and 2D plotting; for 3D plotting, we recommend using the
`old Clawpack MATLAB routines <http://depts.washington.edu/clawpack/users-4.6/matlab_plotting.html>`_
for now.


Basics
=======
Pyclaw includes routines for creating HTML and LaTex plot pages or plotting interactively.
These require a `setplot.py` file that defines the plotting parameters;
see `this help page <http://kingkong.amath.washington.edu/clawpack/users/setplot.html>`_
for more information.  Once you have an appropriate `setplot.py` file,
there are some convenience functions in `$PYCLAW/src/petclaw/plot.py`
for generating these plots.  Assuming you have output files in `./_output`
(which is the default), you can generate HTML plots from Python via

**Comment: Check the plot.py path**

.. doctest::

    >>> from pyclaw import plot
    >>> plot.html_plot() # doctest: +SKIP

This will generate HTML pages with plots and print out a message with the
location of the HTML file.  To launch an interactive plotting session
from within Python, do

.. doctest::

    >>> from pyclaw import plot
    >>> plot.interactive_plot() # doctest: +SKIP

To see a list of commands available in the resulting interactive environment,
type "?".

Plotting result from parallel runs
========================================
By default, when running in parallel, PyClaw outputs data in a binary format.
In order to plot from such files, just replace pyclaw with petclaw in the
commands above; e.g.

.. doctest::

    >>> from petclaw import plot
    >>> plot.interactive_plot() # doctest: +SKIP


More advanced plotting
========================
The plot commands above make use of a file `setplot.py` that is expected to
be found in the current directory and specifies what and how to plot the
results.  If this file is not found, VisClaw will attempt to plot using some
default parameters.  The commands above also assume that the output files
are in a subdirectory of the current directory named `_output`.

For more information on writing your own `setplot.py` file, see the 
`Clawpack documentation <http://depts.washington.edu/clawpack/users-4.6/index.html>`_
and especially the page on `setplot.py <http://depts.washington.edu/clawpack/users-4.6/plotting.html>`_,
as well as the `FAQ <http://depts.washington.edu/clawpack/users-4.6/plotting_faq.html>`_.
