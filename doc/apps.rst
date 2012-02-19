.. _apps:

========================================
Solving other hyperbolic PDEs
========================================
Moving on from the acoustics tutorial, this section explains how to begin
solving other systems of hyperbolic equations with PyClaw.

The built-in PyClaw application scripts
========================================
PyClaw comes with many example application scripts in the directory `pyclaw/apps/`.
These applications are meant to demonstrate the kinds of things that can be done
with PyClaw and are a great place to learn how to use PyClaw.  To run one of them
simply do the following at the command prompt::

    $ cd $PYCLAW/apps/acoustics_1d_homogeneous
    $ python acoustics.py iplot=1

You can run any of the apps similarly by going to the appropriate directory and
executing the Python script.

Command-line options
========================================
For convenience, the scripts are set up to accept certain command-line options.
These usually include the following:

   * use_petsc: set to 1 to run in parallel

   * solver_type: set to classic or sharpclaw

   * iplot: set to 1 to automatically launch interactive plotting after running.
     Note that this shouldn't be used in parallel, as every process will try to plot.

   * htmlplot: set to 1 to automatically creat HTML plot pages after running.

   * outdir: the name of the subdirectory in which to put output files.  Defaults to
     './_output'.

This list of options could easily be extended by modifying the appropriat script.

Adding new applications
========================================
If you have used PyClaw, we'd love to add your application to the built-in scripts.
Please contact us on the `claw-users Google group <http://http://groups.google.com/group/claw-users>`_.
