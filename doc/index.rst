.. _pyclaw:

**********
Quickstart
**********

Install and model a shockwave right away::

    git clone git@github.com:clawpack/clawpack.git
    cd clawpack
    pip install -e .
    cd clawpack/pyclaw/apps/euler_2d
    python shockbubble.py iplot=1

**********
PyClaw is:
**********

    * A **hyperbolic PDE solver** in 1D, 2D, and 3D, including mapped grids and surfaces, built on Clawpack;
    * **Massively parallel** -- the same simple script that runs on your laptop will
      scale efficiently on the world's biggest supercomputers (see :ref:`parallel`);
    * **High order accurate**, with WENO reconstruction and Runge-Kutta time integration
      (see :ref:`solvers`);
    * Simple and intuitive thanks to its Python interface.

PyClaw makes use of the additional Clawpack packages, 
`Riemann <http://github.com/clawpack/riemann>`_ and
`VisClaw <http://github.com/clawpack/visclaw>`_ for Riemann solvers and visualization, 
respectively.

If you have any issues or need help using PyClaw and PetClaw, `send e-mail <claw-users@googlegroups.com>`_
to the `discussion group <http://groups.google.com/group/claw-users>`_.

*************
Documentation
*************

.. toctree::
   :maxdepth: 2

   basics
   going_further
   classes
   for_developers
   troubleshooting
   about
   future


.. _pyclaw_reference:

PyClaw Modules reference documentation
======================================
.. toctree::
   :maxdepth: 1
   
   controller
   solvers
   evolve/limiters
   io
   solution
   state
   geometry
   util

.. _riemann_reference:

Riemann Solvers reference documentation
========================================
The Riemann solvers now comprise a separate package.  For convenience,
documentation of the available pure python Riemann solvers is included
here.  Many other Fortran-based Riemann solvers are available.

.. toctree::
   :maxdepth: 3
   
   rp

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Citing
=======================
If you use PyClaw in work that will be published, please cite the software::

    @misc{pyclaw,
    title={PyClaw software}, 
    url={http://numerics.kaust.edu.sa/pyclaw}, 
    author={Mandli, Kyle T and Ketcheson, David I. and others}, 
    note={Version x.y}
    year={2011}}

and the paper::

    @article{pyclaw-sisc,
            Author = {Ketcheson, David I. and Mandli, Kyle T. and Ahmadia, Aron J. and Alghamdi, Amal and {Quezada de Luna}, Manuel and Parsani, Matteo and Knepley, Matthew G. and Emmett, Matthew},
            Journal = {SIAM Journal on Scientific Computing},
            Month = nov,
            Number = {4},
            Pages = {C210--C231},
            Title = {{PyClaw: Accessible, Extensible, Scalable Tools for Wave Propagation Problems}},
            Volume = {34},
            Year = {2012}}

Please fill in the version number that you used.

If you use the Classic (2nd-order) solver, you may also wish to cite::

    @article{leveque1997,
            Author = {LeVeque, Randall J.},
            Journal = {Journal of Computational Physics},
            Pages = {327--353},
            Title = {{Wave Propagation Algorithms for Multidimensional Hyperbolic Systems}},
            Volume = {131},
            Year = {1997}}

If you use the SharpClaw (high order WENO) solver, you may also wish to cite::

    @article{Ketcheson2011,
            Author = {Ketcheson, D I and Parsani, Matteo and LeVeque, R J},
            Journal = {SIAM Journal on Scientific Computing},
            Number = {1},
            Pages = {A351--A377},
            Title = {{High-order Wave Propagation Algorithms for Hyperbolic Systems}},
            Volume = {35},
            Year = {2013}}


