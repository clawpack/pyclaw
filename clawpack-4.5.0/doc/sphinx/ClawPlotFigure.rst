
.. _ClawPlotFigure:

***************
ClawPlotFigure 
***************


For usage see :ref:`setplot` and :ref:`plotexamples`.

Objects of this class are usually created by the new_plotfigure method of a
:ref:`ClawPlotData` object.

.. class:: ClawPlotFigure


Attributes
==========

  The following attributes can be set by the user:

  .. attribute:: figno : int

    Figure number, by default the next unused number will be used (starting
    at 1001).
    This is usually set as an argument to the new_plotfigure function
    of a :ref:`ClawPlotData` object



  .. attribute:: kwargs : dictionary

    A dictionary of keyword arguments for the figure command.
    For example::

     "{'figsize':[12,5], 'facecolor':[1,0,0]}"

    would specify that the figure should be
    12 inches by 5 inches with a red background.  

    For more options
    see the `matplotlib documentation <http://matplotlib.sourceforge.net/>`_
    for the `figure command
    <http://matplotlib.sourceforge.net/api/figure_api.html#matplotlib.figure.Figure>`_.



  .. attribute:: clf_each_frame : bool

    If True, clear the figure with a clf command at the start of each frame.

  .. attribute:: show : bool

    If False, suppress showing this figure.



Methods
=======

  .. method:: new_plotaxes(name=None)

    Create and return a new object of class :ref:`ClawPlotAxes` associated with this
    ClawPlotFigure object.  A single figure may have several axes on it.

    The name specified is used as a dictionary key.  If None is provided,
    one is generated automatically of the form AXES1, etc.

  .. method:: gethandle()

     Returns the handle for this figure.  

