

.. _python:

***************
Python Hints
***************

.. contents::

Python is a powerful object-oriented interpreted scripting/programming
language. Some version of Python is almost certainly on your computer
already (on unix, linux, OSX type::

    $ python --version 

to find out which version, or just::
    
    $ python

to start a python shell with a prompt that looks like::

    >>>

You may prefer to use `IPython
<http://ipython.scipy.org/moin/>`_, which is a nicer shell
than the pure python shell, with things like command completion and history.
See the `Quick IPython Tutorial
<http://ipython.scipy.org/doc/manual/html/interactive/tutorial.html>`_.

.. _python-install:

Installation of required modules
--------------------------------

To effectively use the pyclaw and Clawpack plotting routines that are
written in Python, you will need version 2.5 or greater
(but **not** 3.0 or above, which is not backwards compatible).  
some modules that are not included in the standard Python distribution. 
Python modules are loaded with the *import* statement in Python and a wide
variety of specialized modules exist for various purposes since people use
Python for many different purposes.

For use with Clawpack, you will need the `Numpy
<http://docs.scipy.org/doc/numpy/user/>`_ module (*Numerical Python*)
that allows working with arrays in much the same way as in Matlab.  
This is distributed as part of 
`SciPy <http://docs.scipy.org/doc/>`_ (*Scientific Python*).

For plotting you will also need the `matplotlib
<http://matplotlib.sourceforge.net/>`_ module which provides Matlab-like
plotting commands for 1d and 2d plots (e.g. contour and pcolor).

Often the easiest way to get all the modules you need is to install the
`Enthought Python Distribution
<http://www.enthought.com/products/epd.php>`_, which is free for academic
users.  Versions are available for Windows, Mac OS X, and Redhat linux.  

For OS X, you might also try the `Scipy Superpack
<http://macinscience.org/?page_id=6>`_.

References and tutorials
------------------------

Some useful links to get started learning Python:

   * `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_
   * `Dive Into Python <http://www.diveintopython.org/>`_

   * `Python tutorial <http://www.python.org/doc/tut/>`_
   * `NumPy User Guide <http://docs.scipy.org/doc/numpy/user/>`_
   * `NumPy for Matlab users <http://www.scipy.org/NumPy_for_Matlab_Users>`_
   * `SciPy Reference Guide <http://docs.scipy.org/doc/scipy/reference/>`_
   * `Matplotlib gallery <http://matplotlib.sourceforge.net/gallery.html>`_
   * `LeVeque's class notes <http://kingkong.amath.washington.edu/uwamath583/sphinx/notes/html/python.html>`_ 
   * `Langtangen's book <http://folk.uio.no/hpl/scripting/>`_ and
     `Introductory slides <http://heim.ifi.uio.no/~hpl/scripting/all-nosplit/>`_


Sage
----

`Sage <http://www.sagemath.org/>`_ is an open source mathematics software 
collection with its own interface and notebook system.  It is based on
Python and the full distribution contains in particular
the Python modules needed for Clawpack.

You can try Sage directly from the
web without downloading, click on "Try Sage Online" link on the
`Sage <http://www.sagemath.org/>`_ webpage.

We are working on incorporating Pyclaw into Sage as an optional package.

