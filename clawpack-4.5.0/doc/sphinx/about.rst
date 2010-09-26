.. _about:

===================
About this software
===================

Clawpack stands for "Conservation Laws Package" and was initially developed
for linear and nonlinear hyperbolic systems of conservation laws, with a
focus on implementing high-resolution Godunov type methods using limiters in
a general framework applicable to many applications.  These finite volume
methods require a "Riemann solver" to resolve the jump discontinuity at the
interface between two grid cells into waves propagating into the neighboring
cells.

Adaptive mesh refinement is included, see :ref:`amrclaw`.

Recent extensions allow the solution of hyperbolic problems that are not in
conservation form.  We are actively working on extensions to parabolic
equations as well.

The "wave propagation" algorithms implemented in Clawpack are discribed in
detail in the book `Finite Volume Methods for Hyperbolic Problems
<http://www.amath.washington.edu/~claw/book.html>`_
Virtually all of the figures in this book were generated using Clawpack and
the source code for each can be found in 
`$CLAW/book <claw/book>`_
See :ref:`book` for a list of available examples with pointers to the codes
and resulting plots.

See the :ref:`biblio` for some pointers to papers describing Clawpack and
the algorithms used in more detail.

A bibliography of older papers using Clawpack can be found 
`here <http://www.amath.washington.edu/~claw/bib.html>`_.  This is out of
date!

.. _license:

License
-------

Clawpack is distributed under the terms of the
Berkeley Software Distribution (BSD) license.  

The licence is in the file `$CLAW/LICENSE.txt <claw/LICENSE.txt>`_ and
reprinted below.

See http://www.opensource.org/licenses/bsd-license.php
for more details.


Copyright (c) 1994--2010, Randall J. LeVeque and others.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 *   Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

 *   Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

 *   Neither the name of the University of Washington nor the names of its
     contributors may be used to endorse or promote products derived from
     this software without specific prior written permission.



THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

.. _authors:

Authors
-------

Many people have contributed to the development of Clawpack since its
inception in 1994.  

Major contributions have been made by the following individuals:

 * Randall J. LeVeque, University of Washington [RJL]

   Originally developed the 1d and 2d Clawpack routines and is the lead
   developer. 

 * Jan Olav Langseth, Norwegian Defence Research Establishment [JOL]

   Developed the 3d algorithms and software with RJL.

 * Marsha Berger, Courant Institute, NYU  [MJB]

   Wrote much of the 2d AMR software originally for the Euler equations and
   worked with RJL to generalize to the AMRCLAW framework.  Later adapted
   this to 3d in work with David McQueen (NYU) and DAC, and continues to be 
   involved in further development of the software.

 * Donna Calhoun, CEA, Paris [DAC]

   Improved Matlab graphics routines and extended to 3d.  Developed 3d AMR
   routines and applications.  Developed ChomboClaw, an extension to allow
   the Clawpack Riemann solvers and frontend to be used together with
   the Chombo AMR package developed by Phil Colella's group at LBL.

 * Sorin Mitran, UNC

   Developed MPI and HDF routines and the Clawpack 4.3 website.  He is also
   the author of BEARCLAW and ASTROBEAR software implementing AMR in a similar
   framework to Clawpack.

 * Peter Blossey, UW

   Developed MPI and HDF versions in Clawpack 4.3.

 * David George, USGS Cascades Volcano Observatory [DLG]

   Developing GeoClaw (with RJL, MJB, KM) for geophysical depth-averaged
   applications.

 * Kyle Mandli, UW

   Main developer of the PyClaw pure python version in directory
   `$CLAW/python/pyclaw <claw/python/pyclaw>`_,
   and the Python/Fortran interface. 

 * Chris Swierczewski, UW

   Working on the Python interface, Sage interface, installation routines.

 * David Ketcheson, KAUST

   Developing SharpClaw for high order methods, e.g. WENO, with Runge-Kutta 
   time stepping.

Numerous students and other users have contributed towards this software, by
finding bugs, suggesting improvements, and exploring its use on new
applications.  Thank you!

.. _funding:

Funding 
-------

Development of this software has been supported in part by

 * NSF Grants DMS-8657319, DMS-9204329, DMS-9303404, DMS-9505021, 
   DMS-96226645, DMS-9803442, DMS-0106511, CMS-0245206,  DMS-0609661,
   DMS-0914942

 * DOE Grants DE-FG06-93ER25181,  DE-FG03-96ER25292, DE-FG02-88ER25053,
   DE-FG02-92ER25139, DE-FG03-00ER2592, DE-FC02-01ER25474

 * AFOSR grant F49620-94-0132, 

 * NIH grant 5R01AR53652-2,

 * ONR grant N00014-09-1-0649

 * The Norwegian Research Council (NFR) through the program no.  101039/420.

 * The Scientific Computing Division at the National Center for Atmospheric
   Research (NCAR).

 * The Boeing Professorship and the Founders Term Professorship in the
   Department of Applied Mathematics, University of Washington.

Any opinions, findings, and conclusions or recommendations expressed in this
material are those of the author(s) and do not necessarily reflect the views
of these agencies. 

.. _citing:

Citing this work
----------------

If you use Clawpack in publications, please cite the following....

   R. J. LeVeque, M. J. Berger, et. al.,  Clawpack Software <version number>,
   `www.clawpack.org <http://www.clawpack.org>`_, <date of access>

Please also cite one of the following regarding the algorithms used in Clawpack:

 * Basic algorithms in 1d and 2d:  [LeVeque97]_, [LeVeque-FVMHP]_

 * 3d algorithms: [LangsethLeVeque00]_

 * AMR: [BergerLeVeque98]_




