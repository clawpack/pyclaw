
.. _ClawPlotData:

**********************
ClawPlotData 
**********************


For usage see :ref:`plotting`, :ref:`setplot`, and :ref:`plotexamples`.


.. class:: ClawPlotData


Attributes
==========


  .. attribute:: outdir : str

    Path to the directory where the Clawpack output is that is to be
    plotted.

  .. attribute:: plotdir : str

    Path to the directory where hardcopy files of plots are to be put.

  .. attribute:: overwrite : bool

    Ok to overwrite old plotdir?  Default is True


  .. attribute:: afterframe : str or function

    Function or string to be executed after producing all plots for each
    frame.  If a string, this string is executed using *exec*.  If a
    function, it should be defined to have a single argument
    "data", [documentation to appear!] For example::

        def afterframe(data):
	    t = data.t
	    print "Done plotting at time %s" % t

  .. attribute:: beforeframe : str or function

    Function or string to be executed before starting to do plots for each
    frame.  If a function, it should be defined to have a single argument
    "data", [documentation to appear!] ::

        def beforeframe(data):


  .. attribute:: printfigs : bool

    True if plots are to be generated and printed (e.g. to png files) before
    making html or latex files.  

    False if the plots have already been generated and the existing versions
    are to be used for making html or latex files.

  .. attribute:: print_format : str
    
    format for hardcopy, default is 'png'


  .. attribute:: print_framenos : list of int's or 'all'

    which frames to print, default is 'all'

  .. attribute:: print_fignos : list of int's or 'all'

    which figures to print for each frame, default is 'all'

  .. attribute:: iplotclaw_fignos : list of int's or 'all'

    which figures to print for each frame in interactive mode, default is 'all'

  .. attribute:: latex : bool

    If True, create a latex file in directory plotdir that to
    display all the plots.

  .. attribute:: latex_fname : str

    Name of latex file, default is 'plots'.  The file created will be e.g.
    plots.tex.

  .. attribute:: latex_title : str

    The title to go at the top of the latex file, default is "Clawpack Results"

  .. attribute:: latex_framesperpage : int or 'all'

    How many frames to try to put on each page.  Default is 'all'.

  .. attribute:: latex_framesperline : int or 'all'

    How many frames to try to put on each line.  Default is 'all'.

  .. attribute:: latex_figsperline : int or 'all'

    How many figures to try to put on each line.  Default is 'all'.
    Recall that several plots may be generated for each frame.

  .. attribute:: latex_pdf : bool

    If True, run pdflatex on the latex file to create e.g., plots.pdf.

  .. attribute:: html : bool

    If True, create a set of html files to display the plots and an index to
    these files called _PlotIndex.html.  These will be in the directory
    specified by the plotdir attribute.



Methods
=======

  .. method:: new_plotfigure(name, figno)

    Create and return a new object of class :ref:`ClawPlotFigure` associated with this
    ClawPlotData object.

  .. method:: getframe(frameno, outdir=None)

    Return an object of class :ref:`ClawSolution` that contains the solution
    read in from the fort.q000N file (where N is frameno).  

    If outdir==None then the outdir attribute of this ClawPlotData
    object is used as the directory to find the data.

    The data, once read in, is stored in a dictionary (the attribute
    framesoln_dict of this ClawPlotData object).  It is read from the fort.q
    file only if it is not already in the dictionary.

    Note that frames read from different outdir's are stored separately (with
    dictionary key (frameno, outdir) if outdir != None).

  .. method:: clearframes(framenos)

    Remove one or more frames from the dictionary framesoln_dict.  
    (Different outdir's not yet implemented.)

  .. method:: clearfigures()

    Clear all plot parameters.  Useful as the first command in a `setplot`
    function to make sure previous parameters are cleared if the file is
    changed and the function is re-executed in an interactive session.


  .. method:: iplotclaw()

    Return True if interactive plotting with iplotclaw is being done.

  .. method:: getfigure(figname)

    Return :ref:`ClawPlotFigure` object with the specified name.

  .. method:: getaxes(axesname, figname=None)

    Return :ref:`ClawPlotAxes` object with the specified name.
    If figname==None then search over all figures and return None if
    it is not found or the name is not unique.

  .. method:: getitem(itemname, axesname=None, figname=None)

    Return :ref:`ClawPlotItem` object with the specified name.
    If axesname==None and/or figname==None 
    then search over all figures and/or axes and return None if
    it is not found or the name is not unique.

  .. method:: showitems()

    Print a list of all the figures, axes, and items defined.

  .. method:: plotframe(frameno)

    Plot all figures for frame number frameno.  Convenience method that
    calls pyclaw.plotters.frametools.plotframe().

  .. method:: printframes()

    Plot and print hardcopy for all frames.  Convenience method that
    calls pyclaw.plotters.frametools.printframes().


.. note:: More methods still to be documented.

