.. _develop:

============================
Information for developers
============================

The PyClaw repository is hosted on Github at 
http://github.com/clawpack/pyclaw.  


Branching development model
============================
PyClaw development follows the model outlined at 
http://nvie.com/posts/a-successful-git-branching-model/, with
one important difference.  The *master* branch is used for development.
New branches are created for releases.
A nice cheat-sheet is available at
http://www.google.com/url?sa=D&q=http://www.globallinuxsecurity.pro/static/git/cheetsheet.pdf.

Contributions in the form of pull requests are also welcome; this approach
is probably the most convenient for occasional contributors, or for major new
contributions like additional solvers.

Bugs
===============
If you find a bug, post an issue with as much explanation as possible on the
Issue tracker at https://github.com/clawpack/pyclaw/issues.  If you're looking 
for something useful to do, try tackling one of the issues listed there.

Developer communication
============================

At the moment, developer communication takes place on the following
google groups:

  * http://groups.google.com/group/pyclaw -- for things relevant only to PyClaw/PetClaw

  * http://groups.google.com/group/claw-dev/ -- for things relevant to the larger Clawpack community

Dependencies
============================

In general, introduction of additional dependencies 
should be avoided.  If you wish to make a change that
will introduce a new dependency (including depending on a more
recent version of a particular package), it should be discussed
on the Google group first.

New versions of existing dependencies will typically be adopted 
either when new functionality provides an important benefit for
PyClaw or when the currently supported version is deemed to be
substantially outdated.

Committing
============================
Always make sure the tests pass before pushing.

Be verbose in your commit messages.

It's helpful to always work in a named branch when
developing a new feature, and to merge with the --no-ff
option so that the history shows distinctly the development
of the feature.


Running the tests
============================
When running the tests, if your machine has multiple cores you can take
advantage of them by doing::

    $ nosetests --processes=2

(replace "2" with the number of processes you want to spawn).

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



