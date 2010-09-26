.. _changes:

##########################
Recent changes
##########################

.. _planned_for_50:

Planned changes for Clawpack 5.0
################################

Substantial changes are planned for Clawpack 5.0.  See the 
`Clawpack wiki <http://kingkong.amath.washington.edu/trac/clawpack>`_ for
details.  Feedback welcome!

You might also view the `tickets <http://kingkong.amath.washington.edu/trac/clawpack/report>`_
on the Trac wiki to see what bugs and features we are working on.

.. _new_in_claw45:

New in Clawpack 4.5
##########################

Overview
--------

For more details, see below.

The svn repository has moved, see below.

The main directory name has changed from `claw` to `clawpack`.

GeoClaw has been added, with a few examples showing how this can
be used for tsunami modeling.  See :ref:`geoclaw_in_45`.

Some new applications have been added.  See :ref:`apps` for a current list
and sample plots.

Subversion repositories and version numbers
-------------------------------------------

We are going to attempt to be more systematic about version numbering
and official releases.  To aid in the long term development of
Clawpack, the Subversion repository at
`http://kingkong.amath.washington.edu/svn/claw4` is being phased out
and replaced by `http://kingkong.amath.washington.edu/svn/clawpack`
The `tags
<http://kingkong.amath.washington.edu/trac/clawpack/browser/tags>`_
subdirectory will contain official releases, starting with the
various versions of Clawpack 4.4 that have been available in the
past as tar files.  These are numbered 4.4.0, 4.4.1, etc.

Note that we have introduced a "micro" version number following the
major and minor version numbers.  Our intention in the future is
to mainly use the micro version number for bug fixes and minor
changes.  New features or more major changes will increment the
minor version number (e.g. going from 4.4 to 4.5).  Major changes to the
structure or functionality will be reflected by incrementing the major
version number.
Contrary to some conventions, we might not always enforce backward
compatibility between minor version numbers.

The `trunk
<http://kingkong.amath.washington.edu/trac/clawpack/browser/trunk>`_ should
be up to date with the most recent release, so that users who want
to use Subversion to keep up to date can check out the trunk and
then use "svn update" to stay current::

    $ svn co http://kingkong.amath.washington.edu/svn/clawpack/trunk  localdir/clawpack

See the `Clawpack wiki`_ for more details.

The `branches <http://kingkong.amath.washington.edu/trac/clawpack/browser/branches>`_
subdirectory contains new development branches, including branches
such as `4.5.x
<http://kingkong.amath.washington.edu/trac/clawpack/browser/branches/4.5.x>`_
for updates that will go into the next release, and branches being
used to develop or test new features.

.. _dir_structure_45:

Directory structure
-------------------

Starting in Version 4.5.0, the main directory is called `clawpack`. 
The location of this directory is where the environment variable
`$CLAW` should point and this convention will be used below.
Within this directory, the structure is currently unchanged.

In Version 5.0 we intend to further rearrange directories.
All the Fortran source code will go in `$CLAW/src` and we plan to introduce
a `$CLAW/lib` for dynamic libraries.  The Makefiles will also change to
reflect these changes.


.. _geoclaw_in_45:

GeoClaw added
-------------

The GeoClaw routines are now incorporated in Clawpack.  
Some documentation is in the section :ref:`geoclaw`.

The main library routines for 2d depth-averaged flow are in
`$CLAW/geoclaw/2d/lib`.

A few examples are in `$CLAW/apps/tsunami`.  See the 
`gallery of sample GeoClaw results <claw/doc/gallery/gallery_geoclaw.html>`_.


The GeoClaw software uses modules and the Makefiles don't always work
properly yet.  If you run into problems, try::

   $ make new

in the applications directory.


.. _new_in_claw44:

New in Clawpack 4.4
##########################

Overview
--------

Clawpack 4.4 consists of the Fortran 77 files from Clawpack 4.3 together
with new Python tools for specifying input data and plotting results.

There is also a preliminary version of a pure Python version of Clawpack,
see :ref:`pyclaw`.


Summary of major changes
------------------------

  * The classic clawpack routines now read data from a file *claw.data*

  * rather than *clawNez.data* and the first line of this file lists the 
    number of space dimensions.   The remainder of the file has the same
    form as before.

  * Rather than modifying *claw.data* it is recommended that you modify
    parameters in the file *setrun.py* and then type 

      $ make .data

    to create the *claw.data* file.  You can modify this to also create
    *setplot.data* or other required data files.

  * The matlab plotting scripts should still work as described in the 4.3
    documentation, but there is now a Python option that uses only open
    source software and provides more powerful plotting tools.
    See :ref:`plotting`.

  * The output routines such as *out1.f*, *out2.f* in classic Clawpack and
    *valout.f* in AMRClaw have been slightly modified to also print ndim to
    the *fort.t* files.  This should not affect Matlab plotting but is
    needed for the new Python plotting routines.

Changes since 4.4.0
-------------------------

 * 10/30/09: Several more changes to Makefiles and amrclaw/2d/lib.

   * This version posted as claw4rev226.tar.gz

   * The subroutines filpatch and prefil are now written as recursive
     subroutines, so filpatch2.f, filpatch3.f, prefil2.f, and prefil3.f have
     been removed.  Also drawrg.f has been removed, an old NCAR graphics
     routine no longer used.  Makefiles in any amr application directory
     will need to be changed to remove these files from the list.

   * For some samples of how to use the latest amrclaw, see e.g.,
     
     * `$CLAW/clawpack/2d/example1/amr  <claw/clawpack/2d/example1/amr/README.html>`_ 
     * `$CLAW/apps/advection/2d/annulus/amr <claw/apps/advection/2d/annulus/amr/README.html>`_ 


 * 10/20/09: Several changes to Makefiles and amrclaw/2d/lib.

   * A new `$CLAW/apps <claw/apps>`_ directory has been added for
     applications.  The ones there now are ones used to debug the amrclaw
     changes, but eventually many more applications from Clawpack 4.3 and
     elsewhere will be put here.

   * New options added to the common Makefile in util/Makefile.common.
     Type "make help" for a list.  Makefiles can now also check for
     dependencies of included files such as call.i used in AMR.

   * New boundary conditions added to amrclaw for problems on the sphere,

   * Dynamic memory for amrclaw - the subroutine init_alloc.f95 was split up into:

    * init_alloc.f90   does initial allocation and the initial size of
      the work array for AMR is specified here.

    * resize_alloc.f90  reallocates for dynamic memory allocation if
      the code runs out of space for AMR.

    * resize_alloc_static.f90  halts with an error message instead of
      reallocating.  For use with compilers that don't support move_alloc,
      such as older versions of gfortran.  This is recommended as the 
      default version in application Makefiles since otherwise it might
      not compile.  Note that some f90 compatible compiler is required
      for using AMR (e.g. gfortran, which is freely available).

    * restart_alloc.f90 is needed when doing a restart with dynamic memory.

    * Note that .f95 files are now relabelled as .f90 since this is
      apparently the standard.

    * Note that Makefiles in user directories that use amrclaw
      will need to be updated to list init_alloc.f90 and
      resize_storage_static.f90 or resize_storage.f90.

   * Several bug fixes in amrclaw/2d/lib


 * 9/18/09: branches/rjl merged back into trunk, includes:

   * Improvements to plotting routines and documentation,

   * More converted examples in the book directory,

   * clawpack/2d/lib directory added with 2d single-grid routines.  

     Similar to version from Clawpack 4.3 but can use setrun.py to set
     runtime parameters and data file is now called claw.data.

     See clawpack/2d/example1 for an example of usage.

   * amrclaw/2d/lib directory added. 
   
     Similar to the version in Clawpack 4.3,
     but with some f95 routines to support dynamic memory allocation.  Also
     gauges are implemented in this version (documentation to appear).

     See clawpack/2d/example1/amr for an example of usage.

