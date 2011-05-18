==========
Plotting
==========

PetClaw output can easily be plotted using the pyclaw plotting routines in
Clawpack.  It's also possible to plot the solution using PETSc.

Plotting with Pyclaw
=====================
Pyclaw includes routines for creating HTML and LaTex plot pages or plotting interactively.
These require a `setplot.py` file that defines the plotting parameters;
see `this help page <http://kingkong.amath.washington.edu/clawpack/users/setplot.html>`
for more information.  Once you have an appropriate `setplot.py` file,
there are some convenience functions in `$PETCLAW/src/petclaw/plot.py`
for generating these plots.  Assuming you have output files in `./_output`
(which is the default), you can generate HTML plots from Python via ::

    >>> from petclaw import plot
    >>> plot.plotHTML()

This will generate HTML pages with plots and print out a message with the
location of the HTML file.  To launch an interactive plotting session
from within Python, do ::

    >>> from petclaw import plot
    >>> plot.plotInteractive()
