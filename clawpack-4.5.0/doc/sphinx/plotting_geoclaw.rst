
.. _plotting_geoclaw:

***************************************
Plotting routines for GeoClaw
***************************************

See :ref:`plotting` for general information about plotting Clawpack results.
GeoClaw results can be viewed using these same tools.  Some addition
functions and useful colormaps are available in the module
`$CLAW/python/pyclaw/plotters/geoplot.py
<claw/python/pyclaw/plotters/geoplot.py>`_

In particular, the following functions are useful to specify as *plot_var*
attributes of a :ref:`ClawPlotItem` (see the module for more documentation):

  topo, land, depth, surface, surface_or_depth

The function *plot_topo_file* is useful for plotting the topography in a
file of the type described in :ref:`topo`.


