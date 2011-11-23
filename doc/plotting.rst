==========
Plotting
==========
Pyclaw includes routines for creating HTML and LaTex plot pages or plotting interactively.
These require a `setplot.py` file that defines the plotting parameters;
see `this help page <http://kingkong.amath.washington.edu/clawpack/users/setplot.html>`_
for more information.  Once you have an appropriate `setplot.py` file,
there are some convenience functions in `$PETCLAW/src/petclaw/plot.py`
for generating these plots.  Assuming you have output files in `./_output`
(which is the default), you can generate HTML plots from Python via ::

    >>> from pyclaw import plot
    >>> plot.html_plot()

This will generate HTML pages with plots and print out a message with the
location of the HTML file.  To launch an interactive plotting session
from within Python, do ::

    >>> from pyclaw import plot
    >>> plot.plot_interactive()

Plotting data files from parallel runs
========================================
By default, when running in parallel, PyClaw outputs data in a binary format.
In order to plot from such files, just replace pyclaw with petclaw in the
commands above; e.g.::

    >>> from petclaw import plot
    >>> plot.plot_interactive()
