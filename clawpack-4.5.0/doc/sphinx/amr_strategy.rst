
.. _amr_strategy:

*****************************************************************
AMR refinement strategy
*****************************************************************

The basic adaptive refinment strategy used in :ref:`amrclaw` is 
to refine on logically rectangular patches.  A single Level 1 grid covers
the entire domain (usually --- if it is too large it may be split into
multiple Level 1 grids).  Some rectangular portions of this grid are covered
by Level 2 grids refined by some refinement factor *R* in each direction
(anisotropic refinement is now allowed too --- see :ref:`setrun_amrclaw`). 
Regions of each Level 2 grid may be covered by Level 3 grids, that are
further refined (perhaps with a different refinement ratio).  And so on.

For the hyperbolic solvers in Clawpack the time step is limited by the
Courant number (see Section :ref:`cfl`), and so if the spatial resolution is
refined by a factor of *R* in each direction then the time step will
generally have to be reduced by a factor *R* as well.  

The AMR code thus proceeds as follows:

 * In each time step on the Level 1 grid(s), the values in all grid cells 
   (including those covered by finer grids) is advanced one time step.
   Before this time step is taken, ghost cells around the boundary of the
   full computational domain are filled based on the boundary conditions
   specified in *$CLAW/amrclaw/Nd/lib/bcNamr.f* (where *N* is the number of
   space dimensions).

 * After a step on the Level 1 grid, *R* time steps must be taken on each
   Level 2 grid, where *R* denotes the desired refinement ratio in
   time from Level 1 to Level 2.  
   
   For each of these time step, ghost cell
   values must be filled in around all boundaries of each Level 2 grid.
   This procedure is defined below in :ref:`amr_bc`.

 * After taking *R* steps on Level 2 grids, values on the Level 1 grid are
   updated to be consistent with the Level 2 grids.  Any cell on Level 1
   that is covered by a Level 2 grid has its *q* value replaced by the
   average of all the Level 2 grid cells lying within this cell.  This gives
   a cell average that should be a better approximation to the true cell
   average than the original value.

 * The updating just described can lead to a change in the total mass
   calculated on the Level 1 grid.  In order to restore global conseravtion,
   it is necessary to do a conservation fix up.  (To be described...)

This style of AMR is often called *Berger-Colella-Oliger* adaptive
refinement, after the papers of Berger and Oliger [BergerOliger84]_ and 
[BergerColella89]_.

The Fortran code in `$CLAW/amrclaw <claw/amrclaw>`_ is based on code
originally written by Marsha Berger for gas dynamics, and merged in Clawpack
in the early days of Clawpack development by MJB and RJL.  The algorithms
used in AMRClaw are described more fully in [BergerLeVeque98]_.


.. _amr_bc:

Ghost cells and boundary conditions for AMR
-------------------------------------------

Consider a Level *k > 1* grid for which we need ghost cells all around the
boundary at the start of each time step on this level.  The same procedure
is used at other levels.

 * Some Level k grids will be adjacent to other Level k grids and so any
   ghost cell that is equivalent to a Level k cell on some other grid has
   values copied from this this grid.

 * Some ghost cells will be in the interior of the full computational domain
   but in regions where there is no adjacent Level k grid.  There will be
   a Level k-1 grid covering that region, however.  In this case the ghost
   cells are obtained by space-time interpolation from values on the Level
   k-1 grid.

 * Some ghost cells will lie outside the full computational domain, where 
   the boundary of the Level k grid lies along the boundary of the full
   domain.  For these cells the subroutine *$CLAW/amrclaw/Nd/lib/bcNamr.f*
   (where *N* is the number of space dimensions) is used to fill ghost cell
   values with the proper user-specified boundary conditions.

For many standard boundary conditions it is not necessary for the user to do
anything beyond setting the *mthbc* parameters in *setrun.py* (see
:ref:`setrun`).  Only if some element of *mthbc* is 0 (indicating
user-specified boundary conditions) is it necessary to modify the library
routine *bcNamr.f* (after copying to your application directory so as not to
damage the library version, and modifying the *Makefile* to point to the new
version).

There some differences between the *bcNamr.f* routine and the *bcN.f*
routine used for the single-grid classic Clawpack routines (which are found in
*$CLAW/clawpack/Nd/lib/bcN.f*).   In particular, it is necessary to check
whether a ghost cell actually lies outside the full computational domain
and only set ghost cell values for those that do.  It should be clear how to
do this from the library version of the routine.

If **periodic boundary
conditions** are specified, this is handled by the AMRClaw software along
with all internal boundaries, rather than in *bcNamr.f*.  With AMR it is not
so easy to apply periodic boundary conditions as it is in the case of a
single grid, since it is necessary to determine whether there is a grid at
the same refinement level at the opposite side of the domain to copy ghost
cell values from, and if so which grid and what index corresponds to the
desired location.  

.. _amr_cluster_fill:

Choosing and initializing finer grids
-------------------------------------

Every few time steps on the coarsest level it is generally necessary to 
revise modify the regions of refinement at all levels, for example to follow
a propagating shock wave.  This is done by

 1. Flagging cells that need refinement according to some criteria.

 2. Clustering the flagged cells into rectangular patches that will form the
    new set of grids at the next higher level.

 3. Creating the new grids and initializing the values of *q* and also any
    *aux* arrays for each new grid.

Clustering is done using and algorithm developed by Berger and Rigoutsis
[BergerRigoutsis91]_ that finds a nonoverlapping set of rectangles that
cover all flagged points and balances the following conflicting goals:

 * Cover as few points as possible that are not flagged,
   to reduce the number of grid cells that must be advanced in each time
   step.

 * Create as few new grids as possible, to minimize the overhead associated
   with filling ghost cells and doing the conservation fix-up around edges
   of grids.

A parameter *cutoff* can be specified (see :ref:`setrun_amrclaw`) to control
clustering.  The algorithm will choose the grids in such a way that at least
this fraction of all the grid points in all the new grids will be in cells
that were flagged as needing refinement.  Usually *cutoff = 0.7* is used, so
at least 70% of all grid cells in a computation are in regions where they
are really needed.

Initializing the new grids at Level k+1 is done as follows:

 * At points where there was already a Level k+1 grid present, this value is 
   copied over.

 * At points where there was not previously a Level k+1 grid, bilinear
   interpolation is performed based on the Level k grids (if the exist at
   this point, if not even coarser grids are used -- **True?**).

.. _amr_flag:

Flagging cells for refinement
-----------------------------

The user can control the criteria used for flagging cells for refinement.

The default procedure is to ...  (explain *tolsp* parameter).

To be continued... describe library routines *allowflag.f* and
*flag2refine.f* and how to modify them.  (But first these should be cleaned
up and regions added to AMRClaw versions!)


.. _regions:

Specifying AMR regions
----------------------

In addition to specifying a tolerance or other criteria for flagging
individual cells as described above, it is possible to specify regions of
the domain so that all points in the region, over some
time interval also specified, will be refined to at least some level
*minlevel* and at most some level *maxlevel*.


**Note:** This is currently available only in :ref:`geoclaw` but should be
carried over to AMRClaw.

This can be automatically specified via parameters set in *setrun.py* (see
:ref:`setrun_regions`). 

To determine whether a grid cell lies in one of the regions specified, the
center of the grid cell is used.  If a mapped grid is being used, the limits
for the regions should be in terms of the computational grid coordinates,
not the physical coordinates.

If a cell center lies in more than one specified region, then the
cell will definitely be flagged for refinement at level k (meaning it should
be covered by a Level k+1 grid) if *k+1 <= minlevel* for any of the regions,
regardless of whether the general flagging criteria hold or not.  
This means the smallest of the various *minlevel* parameters for any region
covering this point will take effect.  Conversely it will definitely **not**
be flagged for refinement if *k+1 > maxlevel* for **all** regions that cover
this point.  This means the largest of the various *maxlevel* parameters for
any region covering this point will take effect.

