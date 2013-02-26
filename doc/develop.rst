.. _develop:

==================================
General information for developers
==================================

The PyClaw repository is hosted on Github at 
http://github.com/clawpack/pyclaw.  


Guidelines for contributing
==================================
When preparing contributions, please follow the guidelines in
:ref:`contribution`.  Also:

    * If the planned changes are substantial or will be backward-incompatible,
      it's best to discuss them on the `claw-dev Google group
      <http://groups.google.com/group/claw-dev>`_ before starting.
      
    * Make sure all tests pass and all the built-in apps run correctly.

    * Be verbose and detailed in your commit messages and your pull request.

    * It may be wise to have one of the maintainers look at your changes before
      they are complete
      (especially if the changes will necessitate modifications of tests
      and/or apps).

    * If your changes are not backward-compatible, your pull request should 
      include instructions for users to update their own application codes.

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

It is almost always okay to add "optional dependencies", for example,  creating
an optional feature or application that is only enabled when the software
dependencies for that feature are satisfied.  :mod:`petclaw` is an example of
this, as it provides an optional means of running many PyClaw applications in
parallel.

In general, introduction of general dependencies 
should be avoided.  If you wish to make a change that
will introduce a new dependency (including depending on a more
recent version of a particular package), it should be discussed
on the `claw-dev Google group`_
first.

New versions of existing dependencies will typically be adopted 
either when new functionality provides an important benefit for
PyClaw or when the currently supported version is deemed to be
substantially outdated.

Running the tests
============================

The PyClaw test suite is built around `nosetests
<http://nose.readthedocs.org/en/latest/>`_ for automatic test discovery, with
supplementary functionality from the :mod:`pyclaw.util` module.  To run the
complete test suite with helpful output, issue the following command at the 
top-level of the pyclaw source directory::

    nosetests -vs

To run the parallel versions of the tests (if petsc4py is installed), run::

    mpirun -n 4 nosetests -vs

(Replace 4 with the number of processes you'd like test on)  

Try prime numbers if you're really trying to break things!

Running serial tests simultaneously
-----------------------------------

When running the tests, if your machine has multiple cores you can take
advantage of them by doing::

    nosetests -vs --processes=2

(replace "2" with the number of processes you want to spawn). However, using
large numbers of processes occasionally causes spurious failure of some tests
due to issues with the operating system.  If you see this behavior, it's best 
to run the tests in serial or with a small number of processes.

Running a specific test
-----------------------

The PyClaw tests are associated with particular applications in the apps/ sub-
directory of the primary repository directory.  If you want to run tests for a
specific application, simply specify the directory containing the application
you are interested in::

   nosetests -vs apps/acoustics_3d_variable

You can also specify a single file to run the tests it contains.

Doctests
--------------

Several of the main PyClaw modules also have doctests (tests in their
docstrings). You can run them by executing the corresponding module::

    cd $PYCLAW/src/pyclaw
    python grid.py
    python state.py

If the tests pass, you will see no output.  You can get more output by using 
the `-v` option::

    python state.py -v

Writing New Tests
==================

If you contribute new functionality to PyClaw, it is expected that you will also
write at least one or two new tests that exercise your contribution, so that
further changes to other parts of PyClaw or your code don't break your feature.
You do not have to use any of the functionality offered by pyclaw.util, but it
may simplify your test-writing and allow you to check more cases than you would
easily specify by hand.

The most important function in :mod:`pyclaw.util` is
:func:`pyclaw.util.gen_variants`, which allows you to perform combinatorial
testing without manually specifying every feature you'd like to perform.
Currently, :func:`~pyclaw.util.gen_variants` can multiplicatively exercise
kernel_languages (Fortran or Python) and pure PyClaw or PetClaw implementations.
This allows you to write one function that tests four variants.

Another function provided by :mod:`~pyclaw.util` is
:func:`pyclaw.util.test_app`. The :func:`~pyclaw.util.test_app` function will
run an application as if started from the command line with the specified
keyword arguments passed in.  This is useful for testing specific code that does
not necessarily work with :mod:`petclaw`, for example, and is not expected to.

You will notice that both :func:`~pyclaw.util.gen_variants` and
:func:`~pyclaw.util.test_app` require a `verifier` method as an argument. 
These functions both effectively run tests and verify output with the following
function calls::
 
        output = application(**kwargs)
        check_values = verifier(output)

The `verifier` method needs to return `None` if there is no problem with the
output, or a sequence of three values describing what was expected, what it
received, and more details about the error.  A very simple `verifier` method
that you can use is :func:`pyclaw.util.check_diff`, which can use either an
absolute tolerance or a relative tolerance to compare an expected value against
the test output from the application.

See apps/acoustics_1d_homogeneous/test_acoustics.py for a comprehensive example
of how to use :func:`~pyclaw.util.gen_variants` and
:func:`~pyclaw.util.check_diff`. See apps/shallow_sphere/test_shallow_sphere.py
for an example that uses :func:`~pyclaw.util.test_app` and also loads a known
solution from disk using numpy.

Catching errors with Pyflakes and Pylint
===========================================

Pyflakes and Pylint are Python packages designed to help you catch errors or
poor coding practices.  To run pylint on the whole PyClaw package, do::

    cd $PYCLAW
    pylint -d C pyclaw

The `-d` option suppresses a lot of style warnings, since PyClaw doesn't
generally conform to PEP8.  To run pylint on just one module, use something
like::

    pylint -d C pyclaw.state

Since pylint output can be long, it's helpful to write it to an html file
and open that in a web browser::

    pylint -d C pyclaw.state -f html > pylint.html

Pyflakes is similar to pylint but aims only to catch errors.  If you
use Vim, there is a nice extension package 
`pyflakes.vim <https://github.com/kevinw/pyflakes-vim>`_
that will catch errors as you code and underline them in red.

Checking test coverage
========================
You can use nose to see how much of the code is covered by the current
suite of tests and track progress if you add more tests ::

    nosetests --with-coverage --cover-package=pyclaw --cover-html

This creates a set of html files in `./cover`, showing exactly which lines
of code have been tested.
