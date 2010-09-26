

.. _plotting_faq:

***********************
Plotting hints and FAQ
***********************

.. seealso:: 
    :ref:`plotting`, 
    :ref:`setplot`, 
    :ref:`plotexamples` 


.. contents::

More to come!


How to plot a something other than a component of q?
----------------------------------------------------

Objects of class :ref:`ClawPlotItem` have an attribute ``plot_var``.  If
this is set to an integer than the corresponding component of q is plotted
(remembering that Python indexing starts at 0, so ``plot_var = 0`` will
specify plotting the first component of q, for example).

If plot_var is a function then this function is applied to applied to 
:ref:`current_data` and should return the array of values to be plotted.
For an example, see :ref:`plotexample-acou-1d-6`.

Sometimes you want to plot something other than the solution on the grid, 
for example to add another feature to a plot of the solution. This can be
done via an ``afteraxes`` command, for example, which is called after all
items have been plotted on the current axes.  See :ref:`ClawPlotAxes` for
details and an example.

How to specify ``outdir``, the directory to find ``fort.*`` files for plotting?
-------------------------------------------------------------------------------

This is normally determined by the ``outdir`` attribute of
the :ref:`ClawPlotData` object directing the plotting.  But see the next FAQ
for the option of using different directories for some plot items (e.g. to
compare results of two computations).

If you are making a set of hardcopy plots using::

    $ make .plots

then ``outdir`` is specified in the Makefile by setting the ``CLAW_OUTDIR``
variable.

If you are making plots interactively using Iplotclaw_, then you can
directly specify the ``outdir`` as a parameter, e.g.::

    In[1]: ip=Iplotclaw(outdir="_output");   ip.plotloop()

If you don't specify this parameter, `Iplotclaw`_ will look for a file
``.output`` in the current directory.  If you created the ``fort.*`` files by
the command::

    $ make .output

then the output directory is set in the Makefile and the file ``.output``
contains the path to the output directory.

If the file ``.output`` does not exist,  ``outdir = '.'`` is used by
default, the current directory.  

Note that if you stop a calculation mid-stream using ``<ctrl>-C``, the file
``.output`` may not exist or be correct, since this file is written after
the execution finishes.  

How to specify a different ``outdir`` for some plot item?
-------------------------------------------------------------

If you want one plot item on an axis to use the default ``plotdata.outdir``
while another to take data from a different directory (in order to compare
two computations, for example), you can set the ``outdir``
attribute of a :ref:`ClawPlotItem` directly.  If you do not set it, by
default it inherits from the :ref:`ClawPlotFigure` object this item belongs
to.

For example, you might have the following in your ``setplot`` function::

    plotfigure = plotdata.new_plotfigure(name='compare', figno=1)
    plotaxes = plotfigure.new_plotaxes()

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    import os
    plotitem.outdir = os.path.join(os.getcwd(), '_output2')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-+'
    plotitem.color = 'r'

This would plot results from ``plotdata.outdir`` as blue circles and results
from ``./_output2`` as red plus signs.  It's best to give the full path
name, e.g. as done here using ``os.path.join(os.getcwd(), '_output2')``.



How to debug setplot.py?
--------------------------

Suppose you are working in an interactive Python shell such as ipython and
encounter the following when trying to plot with `Iplotclaw`_::

    In [3]: ip=Iplotclaw(); ip.plotloop()
    *** Error in call_setplot: Problem executing function setplot
    *** Problem executing setplot in Iplotclaw
        setplot =  setplot.py
    *** Either this file does not exist or 
        there is a problem executing the function setplot in this file.
    *** PLOT PARAMETERS MAY NOT BE SET! ***
    
    Interactive plotting for Clawpack output... 
    
    Plotting data from outdir =  _output
    Type ? at PLOTCLAW prompt for list of commands
    
        Start at which frame [default=0] ? 
    
    
This tells you that there was some problem importing ``setplot.py``, but is not
very informative and it is hard to debug from within the
``Iplotclaw.plotloop``
method. You may also run into this if you modify ``setplot.py``
(inadvertantly introducing a bug)
and then use the ``resetplot`` option::

    PLOTCLAW > resetplot
    Executing setplot from  setplot.py
    *** Error in call_setplot: Problem executing function setplot
    *** Problem re-executing setplot
    PLOTCLAW > 


If you can't spot the bug by examing ``setplot.py``, it is easiest to debug
by exiting the plotloop and doing::
    
    PLOTCLAW > q
    quitting...
    
    In [4]: import setplot
    In [5]: pd = ip.plotdata
    In [6]: pd = setplot.setplot(pd)
    ---------------------------------------------------------------------------
    AttributeError                            Traceback (most recent call last)
    
          8 
          9     # Figure for q[0]
    ---> 10     plotfigure = plotdata.new_plotfgure(name='q[0]', figno=1)
         11 
         12     # Set up for axes in this figure:
    
    AttributeError: 'ClawPlotData' object has no attribute 'new_plotfgure'
    
    
In this case, the error is that ``new_plotfigure`` is mis-spelled.

In ipython you can also easily turn on the Python debugger pdb::

    In [9]: pdb
    Automatic pdb calling has been turned ON

    In [10]: pd = setplot.setplot(pd)
    ---------------------------------------------------------------------------
    AttributeError                            Traceback (most recent call last)
          8 
          9     # Figure for q[0]
    ---> 10     plotfigure = plotdata.new_plotfgure(name='q[0]', figno=1)
         11 
         12     # Set up for axes in this figure:

    AttributeError: 'ClawPlotData' object has no attribute 'new_plotfgure'

    ipdb> 

For more complicated debugging you could now explore the current state using
any pdb commands, described in the `documentation
<http://docs.python.org/library/pdb.html>`_.  See also 
the `ipython documentation
<http://ipython.scipy.org/doc/manual/html/index.html>`_.
