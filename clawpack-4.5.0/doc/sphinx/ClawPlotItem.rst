.. _ClawPlotItem:

************
ClawPlotItem
************


For usage see :ref:`setplot` and :ref:`plotexamples`.

Objects of this class are usually created by the new_plotitem method of a
:ref:`ClawPlotAxes` object.

The examples shown below are fragments of a typical *setplot*
function which assume *plotaxes* is an instance of :ref:`ClawPlotAxes`.

.. seealso:: :ref:`setplot` and :ref:`plotexamples`
.. seealso:: :ref:`plotexamples`

.. class:: ClawPlotItem

Attributes
==========

The following attributes can be set by the user:

.. attribute:: plot_type : str

    Type of plot desired, one of the following:

    * '1d_plot' : one dimensional line or set of points plotted using the
      matplotlib plot command.
    * '1d_from_2d_data' : 1d plot generated from 2d data, for example as a
      slice of the data or a scatter plot of data that should be radially
      symmetric,
    * '1d_fill_between' : 1d filled plot between two variable specified by
      the attributes *plot_var* and *fill_var2*.

    * '2d_contour' : two dimensional contour plot,
    * '2d_pcolor' : two dimensional pcolor plot,
    * '2d_schlieren' : two dimensional Schlieren plot,
    * '2d_grid' : two dimensional plot of only the grids, no data

.. attribute:: outdir : str or None

     Directory to find data to be plotted.  
     If None, then data comes from the outdir attribute of the parent 
     ClawPlotData item.

.. attribute:: plot_var : int or function

     If an integer, then this specifies which component of q to plot (with
     the Python convention that plot_var=0 corresponds to the first
     component).

     If a function, then this function is applied to q on each grid to
     compute the variable var that is plotted.  The signature is

     * var = plot_var(current_data)

     The :ref:`current_data` object holds data from the current frame that can be
     used to compute the variable to be plotted.  


.. attribute:: afteritem : str or function or None

     A string or function that is to be executed after plotting this item.
     If a string, this string is executed using *exec*.  If a
     function, it should be defined to have a single argument
     :ref:`current_data`.
     
     For example::

        def afteritem(current_data):

.. attribute:: aftergrid : str or function or None

     A string or function that is to be executed after plotting this item on
     each grid. (There may be more than 1 grid in an AMR calculation.)
     If a string, this string is executed using *exec*.  If a
     function, it should be defined to have a single argument
     "data", [documentation to appear!] 
     
     For example::

        def aftergrid(current_data):
            cd = current_data
	    print "On grid number %s, xlower = %s, ylower = %s" \
	          % (cd.gridno, cd.xlower, cd.ylower)
     
     would print out the grid number and lower left corner for each grid in
     a 2d computation after the grid is plotted.



.. attribute:: MappedGrid : bool

     If True, the grid mapping specified by the *mapc2p* attribute of the
     underlying `ClawPlotData` object should be applied to the grid before
     plotting.


.. attribute:: show : bool

     If False, plotting of this object is suppressed.



**The other attributes required depend on the plot_type, as summarized
below:**

Special attributes for all 1d plots,  plot_type = '1d...'
----------------------------------------------------------

.. attribute:: plotstyle : str

     Anything that is valid as a fmt
     group in the `matplotlib plot command
     <http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot>`_.
     For example:

     * '-' for a solid line, '- -' for a dashed line,
     * 'o' for circles, 'x' for x's, '-o' for circles and a line,
     * 'bo' for blue circles (though if the *color* attribute is also set
       that will overrule the color in the format string).

.. attribute:: color : str

     Any matplotlib color, for example red can be specified as 'r' or 'red'
     or '[1,0,0]' or '#ff0000'.

Special attributes for plot_type = '1d_plot'
----------------------------------------------------------

No extra attributes.

Special attributes for plot_type = '1d_fill_between'
----------------------------------------------------------

This gives a filled polygon between two curves using the `matplotlib
fill_between command <http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.fill_between>`_.

.. attribute:: plot_var : int or function

    as described above defines one curve,


.. attribute:: plot_var2 : int or function

    defines the second curve for the fill_between command. 
    Default is the zero function.

.. attribute:: fill_where : str or None

    defines the *where* attribute of the fill_between command.

  Example::

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = 0    # means use q[:,0] 

  would produce a filled curve between y=q[:,0] and y=0.  

  Example::

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = 0    # means use q[:,0] 
    plotitem.plot_var2 = 1

  would produce a filled curve between y=q[:,0] and y=q[:,1].
    
.. _1d_from_2d_data:

Special attributes for plot_type = '1d_from_2d_data'
------------------------------------------------------

.. attribute:: map_2d_to_1d : function

  Example:  In a 2d computation where the solution q[:,:,0] should be
  radially symmetric about (x,y)=(0,0), the following will result in a
  scatter plot of the cell values q[i,j,0] vs. the radius r(i,j)::

    def q0_vs_radius(current_data):
        # convert 2d (x,y,q) into (r,q) for scatter plot
        from numpy import sqrt
        x = current_data.x
        y = current_data.y
        r = sqrt(x**2 + y**2)
        q0 = current_data.var   # the variable specified by plot_var
        # q0 = current_data.q[:,:,0]   # would also work
        return r,q0

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = 0     # use q[:,:,0]
    plotitem.plotstyle = 'o'  # symbol not line is best for scatter plot
    plotitem.map_2d_to_1d = q0_vs_radius   # the function defined above
    
  See :ref:`current_data` for a description of the *current_data* argument.


Special attributes for all 2d plots,  plot_type = '2d...'
------------------------------------------------------------

.. attribute:: gridlines_show : bool

     If True, draw the grid lines on the plot.  
     The attribute 'amr_gridlines_show' should be used for AMR computations
     to specify that gridlines should be shown on some levels and not
     others. See :ref:`amr_attributes`.

.. attribute:: gridedges_show : bool

     If True, draw the edges of grids, mostly useful in AMR computations.

Special attributes for plot_type = '2d_contour'
------------------------------------------------------

.. attribute:: contour_levels : numpy array or None

     If a numpy array, the contour levels.  If None, then the next three
     attributes are used to set the levels.

.. attribute:: contour_nlevels : int

     Number of contour levels

.. attribute:: contour_min : float

     Minimum contour level

.. attribute:: contour_max : float

     Maximum contour level

.. attribute:: contour_colors : color specification

     Colors of contour lines.  Can be a single color such as 'b' or
     '#0000ff', or a colormap.

.. attribute:: amr_contour_colors : list of color specifications

     As with other attributes (see :ref:`amr_attributes` below), 
     instead of contour_colors you can specify
     *amr_contour_colors*
     to be a list of colors (or colormaps) to use on each AMR level, e.g.::
       
         amr_contour_colors = ['k','b','r']

     to use black lines on Level 1, blue on Level 2, and red for all
     subsequent levels.  This is useful since with the matplotlib contour
     plotter you will see both fine and course grid lines on top of one
     another in refined regions (Matplotlib lacks the required
     hidden line removal to blank out the lines from coarser grids easily.
     See also the next attributes.)

.. attribute:: contour_show : boolean

     Show the contour lines only if this attribute is true.  This is most
     commonly used in the form of the next attribute,


.. attribute:: amr_contour_show : list or tuple of booleans

     Determines whether to show the contour lines on each AMR level.  Useful
     if you only want to view the lines on the finest grids.


.. attribute:: contour_kwargs : dictionary

     Other keyword arguments for the contour command.

Special attributes for plot_type = '2d_pcolor'
-------------------------------------------------

.. attribute:: pcolor_cmap : matplotlib colormap

.. attribute:: pcolor_cmin : float

.. attribute:: pcolor_cmax : float

     In general you should specify *pcolor_cmin* and *pcolor_cmax* to
     specify the range of q values over which the colormap applies.  If they 
     are not specified they will be chosen automatically and may vary from
     frame to frame.  Also, if AMR is used, they may vary from grid to grid,
     yielding very confusing plots.

.. attribute:: pcolor_colorbar : bool

     If True, a colorbar is added to the plot.


.. _amr_attributes:

AMR Attributes
==============

Many attributes listed above also have a second related attribute with the
same name pre-pended with *amr_*.  If this attribute is set, it should be a
list whose elements are of the type specified by the original name, and the
elements of the list will be used for different AMR refinement levels.

For example, the following commands::

    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.contour_color = 'r'

will result in all contour lines being red on all levels of AMR.  On the
other hand::

    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.amr_contour_color = ['k', 'b']

will result in contour lines on grids at level 1 being black and on grids of
level 2 or higher being blue.  

Note that if the list is shorter than the number of levels, the last element
is used repeatedly.

If both attributes *contour_color* and *amr_contour_color* are set, 
only *amr_contour_color* is used.

A common use is to show grid lines only on coarse levels, not on finer
levels, e.g.::

    plotitem.amr_gridlines_show = [1,1,0]

will result in gridlines being shown only on levels 1 and 2, not on finer
levels.




Methods
=======

.. method:: getframe(frameno) 

     Returns an object of class :ref:`Solution` 
     containing the solution being
     plotted by this object for frame number frameno.

.. method:: gethandle()

     Returns the handle for this item.  

