.. _develop:

============================
Information for developers
============================

The PyClaw repository is hosted on Github at 
http://github.com/clawpack/pyclaw.  


Guidelines for contributing
==================================
When preparing contributions, please follow the guidelines in
:ref:`using-git`.  Also:

    * If the planned changes are substantial or will be backward-incompatible,
      it's best to discuss them on the `claw-dev Google group
      <http://groups.google.com/group/claw-dev`_ before starting.
      
    * Make sure all tests pass and all the built-in apps run correctly.

    * Be verbose and detailed in your commit messages and your pull request.

    * It may be wise to have one of the maintainers look at your changes before
      they are complete
      (especially if the changes will necessitate modifications of tests and/or apps).

    * If your changes are not backward-compatible, your pull request should include
      instructions for users to update their own application codes.

Bugs
===============
If you find a bug, post an issue with as much explanation as possible on the
Issue tracker at https://github.com/clawpack/pyclaw/issues.  If you're looking 
for something useful to do, try tackling one of the issues listed there.

Developer communication
============================

As PyClaw is part of the family of Clawpack codes, developer communication
takes place on the google group at http://groups.google.com/group/claw-dev/.

Dependencies
============================

In general, introduction of additional dependencies 
should be avoided.  If you wish to make a change that
will introduce a new dependency (including depending on a more
recent version of a particular package), it should be discussed
on the `claw-dev Google group <http://http://groups.google.com/group/claw-dev>`_.
first.

New versions of existing dependencies will typically be adopted 
either when new functionality provides an important benefit for
PyClaw or when the currently supported version is deemed to be
substantially outdated.


Running the tests
============================
When running the tests, if your machine has multiple cores you can take
advantage of them by doing::

    $ nosetests --processes=2

(replace "2" with the number of processes you want to spawn).
However, using large numbers of processes occasionally causes spurious failure
of some tests due to issues with the operating system.  If you see this
behavior, it's best to run the tests in serial or with a small number of
processes.

It is also possible to perform only a subset of the regression tests
(e.g. pure python code or python and fortran code, classic clawpack or
sharpclaw solver, explicit or implicit time stepping, etc.). This can be
accomplished by passing some attributes to nose. The attributes are already
defined in the regression tests suite and they are:

    * solver_type: classic or sharpclaw
    * kernel_language: python or fortran
    * petsc: True or False
    * time_stepping_mode: explicit
    * time_stepping_method: ForwardEuler, SSP33, SSP104 
    * speed: fast or slow

The attribute 'time_stepping_method' is only used in combination with
'solver_type = sharpclaw' because the classic clawpack implements the
Lax-Wendroff scheme.

The attributes can be used in the following ways:

    * Logic AND: run only the regression tests that have the listed attributes ::
    
        $ nosetests -a attribute-1 = value-1,attribute-2 = value-2,attribute-3 = value-3

    * Logic OR: run the regression tests that have at least one of the listed attributes :: 
    
        $ nosetests -a attribute-1 = value-1 -a attribute-2 = value-2 -a attribute-3 = value-3



