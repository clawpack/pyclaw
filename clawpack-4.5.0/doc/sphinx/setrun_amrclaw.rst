

.. _setrun_amrclaw:

*****************************************************************
Specifying AMRClaw parameters in `setrun.py`
*****************************************************************

Since AMRClaw is an extension of Clawpack, all of the parameters that
are required for Clawpack are also needed by AMRClaw.  See
:ref:`setrun` for a discussion of these and a 
description of `setrun.py` input scripts more generally.

In addition, a number of other parameters should be set in the `setrun.py`
file in any AMRClaw application.

It is best to look at a specific example while reading this section, for
example  :ref:`setrun_amrclaw_sample`.

The function `setrun` in this module is essentially the same as for Clawpack,
except that it expects to be called with *claw_pkg = 'amrclaw'*.  This call
should be performed properly by the Makefile if you have *CLAW_PKG =
amrclaw* set properly there.

The new parameter in this module start at 
:ref:`setrun_amrclaw_sample_parameters` in the sample file.

A brief summary of these:

.. attribute:: mxnest : integer

   **mxnest** is the maximum number of refinement levels to use.  
   *mxnest=1* corresponds to a single grid run and should give essentially the
   same results as the classic Clawpack would give (not identical probably
   because different routines are used with minor variations in
   implementation). Checking that this works is a good first step
   in converting a code to \amrclaw.

   *mxnest > 1* then more then one level is used.  

   *mxnest < 0* means *abs(mxnest)* levels are used.  The negative sign
   indicates that anisotropic refinement may be used, which affects the next
   parameters.

   This will be cleaned up in Clawpack 5.0.

.. attribute:: inratx : list of integers

   *inratx* is a list of refinement ratios in the x direction.  
   *inratx[k] = R* means that in refining from level k+1 to k+2 the x
   direction will be refined by a factor R.  (Here Python indexing starting
   at 0 is used, so *inratx[0]* is the ratio from Level 1 to Level 2.)

   If *mxnst > 0* then *inratx* also determines the refinement ratio in y
   and t (and in z for 3d problems).  

.. attribute:: inraty, inratz, inratt : lists of integers

   *inraty*, *inratz*, *inratt* are only used if *mxnest < 0* and are lists
   of refinement ratios in the respective directions.

.. attribute:: auxtype : list of strings

   If *maux > 0* then for each component of *aux* there should be a
   corresponding element of *auxtype* from the list below:

     'xleft'
        a value associated with the left edge of a cell in the x-direction

     'yleft'
        a value associated with the left edge of a cell in the y-direction

     'zleft'
        a value associated with the left edge of a cell in the z-direction

     'center'
        a value associated with a cell center

     'capacity'
        a cell-centered capacity function


   The *auxtype* array is required for adaptive refinement because
   auxiliary arrays must be handled slightly differently at refinement
   boundaries depending on how these values are used.

   A cell-centered auxiliary value such as the density or impedance in a 
   heterogeneous acoustics problem would have type *center*. 
   On the other hand, in a variable-coefficient advection problem we may
   want to store the normal velocity at each edge of the cell. In two
   dimensions we might use  one component of *aux* to store a left-edge
   velocity in the x-direction and another to store the left-edge (i.e.
   bottom) velocity in the y-direction.

   At most one component may have type *'capacity'*, and the value of
   *mcapa* described in :ref:`setrun` should be set in a consistent
   manner (modulo the fact that Fortran indexing starts at 1 and Python at
   0, so if *mcapa = k* then *auxtype[k-1] = 'capacity'* should be set).

.. attribute:: checkpt_iousr : integer

   A  A checkpoint file is dumped every *checkpt_iousr* time
   steps on the coarse grid.  These are binary files with names of the form
   *fort.chkXXXX* where *XXXX* is the step number.   

   **Note:** these files are typically very large!

   The solution and grid structure is printed out in a form that can be used
   to later restart the calculation  from this point.   This is useful when
   doing long runs in case the computer goes down or the algorithm fails at some
   point in the calculation.  It is also useful if you want to go to some
   large time and then start doing frequent outputs in order to examine the
   time-evolution of the solution more carefully.

   In addition to creating a checkpoint file every *checkpt_iousr* time steps, a
   final checkpoint file is created at the end of the computation.  This can be
   used to restart the calculation from the final time if you wish to evolve it
   further.  Setting *checkpt_iousr* to a sufficiently
   large integer will cause a checkpoint
   file to be written only at the end of the computation.

   If *checkpt_iousr = 0* then no checkpoint files are generated, not even at
   the end.

   If *checkpt_iousr < 0* then the attribute *tchk* should also be set, to a 
   list of *abs(checkpt_iousr)* times  when checkpoint files are desired.

   This will be cleaned up in Clawpack 5.0.

.. attribute:: restart : boolean

   If *restart = True* then a restart is performed.
   Information read in from the file *restart.data* is used to resume a
   previous calculation.  An appropriate checkpoint file *fort.chkXXXX*
   should be renamed *restart.data* in order to use this option.
   (And generally moved from the *_output* directory to the directory from
   which the code is being run.)

   When a restart is performed, other parameters in this *amr2ez.data* file
   should be consistent with values used in the previous calculation, with
   some exceptions:

     * The final time *tfinal* can be increased,
     * Others?

   Note that when restarting, the output files will continue to be numbered
   consecutively from the previous run.

   
.. attribute:: tol: float

   Error tolerance used in Richardson error estimation.  Cells are
   flagged for refinement if the error estimate is greater than *tol*.
   Richardson estimation requires taking two time steps on the current grid
   and comparing the result with what's obtained by taking one step on a
   coarsened grid.  

   If *tol < 0*, Richardson estimation is not used.  

.. attribute:: tolsp: float

   Error tolerance used in a simpler approach of estimating the spatial
   gradient of the solution and flagging points where this estimate is
   larger than *tolsp*.  See :ref:`amr_strategy` for more information.

.. attribute:: kcheck : int

   How often to regrid: error estimation and regridding is performed every
   *kcheck* time steps on each level.
   
.. attribute:: ibuff : int

   Size of the buffer zone around flagged cells.
   Certain cells are flagged for refinement and then clustered (see
   :ref:`amr_strategy`) into finer grids.  In addition to the cells flagged by the
   error estimation, all cells within *ibuff* cells of these are also
   flagged.  This insures that structures in the solution that require
   refinement will remain in the refined region for at least *ibuff* time
   steps, since the Courant number must be no greater than 1.  The value of
   *ibuff* should generally be consistent with the value of *kcheck*,
   with *ibuff >= kcheck* if the Courant number is close to 1.
   
.. attribute:: cutoff : float

   Parameter used in the clustering algorithm (see
   :ref:`amr_strategy`).  Typically 0.7 is a good value.


