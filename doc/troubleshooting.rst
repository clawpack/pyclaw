.. _troubleshooting:

********************
Troubleshooting
********************

This page lists some of the most common difficulties encountered in 
installing and running PyClaw.  If you do not find a solution for your
problem here, please e-mail the 
`claw-users Google group <http://http://groups.google.com/group/claw-users>`_.
You may also wish to consult the `list of known issues <https://github.com/clawpack/pyclaw/issues>`_.

Compilation errors
********************
Two frequent sources of compilation error are:

    * Your environment variable FC is set to g77 or another Fortran 77 compiler.
      FC should be undefined or set to a Fortran 90 compiler.
      If you have installed gfortran, you could set::

        $ export FC=gfortran

     in your .bash_profile.

    * Conflicts between 32-bit and 64-bit files.  This has been encountered on
      Mac OS X with 32-bit Enthought Python.  We recommend using a 64-bit Python
      install, such as that available from Enthought (free for academics).
      The 32-bit EPD has also been known to cause a plotting issue with PyClaw
      in which plotting becomes extremely slow.
