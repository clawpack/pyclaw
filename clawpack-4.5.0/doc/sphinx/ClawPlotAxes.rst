
.. _ClawPlotAxes:

**************************************
ClawPlotAxes 
**************************************


For usage see :ref:`setplot` and :ref:`plotexamples`.

Objects of this class are usually created by the new_plotaxes method of a
:ref:`ClawPlotFigure` object.

.. seealso:: :ref:`setplot` and :ref:`plotexamples`

.. class:: ClawPlotAxes


Attributes
==========

  The following attributes can be set by the user:

  .. attribute:: name : str

  .. attribute:: title : str

    The title to appear at the top of the axis.  Defaults to string
    specified by the *name* attribute.

  .. attribute:: title_with_t : bool

    If True, creates a title of the form "%s at time t = %s" % (title, t)

  .. attribute:: axescmd : str

  The command to be used to create this axes, for example:
    *  "subplot(1,1,1)" for a single axes filling the figure
    *  "subplot(2,1,1)" for the top half
    *  "axes([0.1, 0.1, 0.2, 0.8])" for a tall skinny axis.

  See the matplotlib documentation for axes.

  .. attribute:: xlimits : array [xmin, xmax]  or 'auto'

     The x-axis limits if an array with two elements, or choose
     automatically

  .. attribute:: ylimits : array [ymin, ymax]  or 'auto'

     The y-axis limits if an array with two elements, or choose
     automatically

  .. attribute:: afteraxes : function or str

     A string or function that is to be executed after creating all 
     plot items on this axes.
     If a string, this string is executed using *exec*.  If a
     function, it should be defined to have a single argument
     :ref:`current_data`.  

     The string version is useful for 1-liners such as::

        afteraxes = "pylab.title('My custom title')"

     pylab commands can be used, since pylab has been imported into the
     plotting module.
     
     The function form is better if you want to do several things, or if you
     need access to the data stored in :ref:`current_data`.  For example::

        def afteraxes(current_data):
            # add x- and y-axes to a 1d plot already created
            from pylab import plot

            xlower = current_data.xlower
            xupper = current_data.xupper
            plot([xlower, xupper], [0.,0.], 'k')   # x-axis

            # Get y limits from variable just plotted, which is
            # available in current_data.var.  
            ymin = current_data.var.min() 
            ymax = current_data.var.max()
            plot([0.,0.], [ymin,ymax], 'k')  # y-axis


  .. attribute:: show : bool

     If False, suppress showing this axes and all items on it.


Methods
=======

  .. method:: new_plotitem(name=None, plot_type)

     Returns an object of class :ref:`ClawPlotItem` associated with this axes.
     A single axes may have several items associated with it.

     The name specified is used as a dictionary key.  If None is provided, 
     one is generated automatically of the form ITEM1, etc.


  .. method:: gethandle()

     Returns the handle for this axes.  

