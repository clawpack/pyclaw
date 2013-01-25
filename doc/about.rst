.. _about:


=======================
About PyClaw
=======================
PyClaw is an open-source, free project.  If you find PyClaw useful,
please let us know via e-mail to pyclaw@googlegroups.com.

Contributors
=======================
Many people have contributed to PyClaw, some of them very substantial parts of
the package. Their work is greatly appreciated: no open source project can
survive without a community. The following people contributed major parts of
the library (in alphabetical order)

    * Aron Ahmadia: PETSc integration; I/O; programmatic testing framework; general design.

    * Amal Alghamdi: Initial development of PetClaw, many bug-fixes and enhancements.

    * Jed Brown: Implicit time stepping (still experimental).

    * Lisandro Dalcin: Fortran wrapping; PETSc integration; general efficiency.

    * Matthew Emmett: PyWENO integration.

    * David Ketcheson: General maintenance and development; incorporation of SharpClaw routines.

    * Matthew Knepley: General design; PETSc integration.

    * Grady Lemoine: Interleaving and cache-optimization of 3D Classic routines.

    * Kyle Mandli: Initial design and implementation of the PyClaw framework.

    * Matteo Parsani: Mapped grids; Python-Fortran interfacing; implicit time stepping.

Further contributions to the package are most welcome.  If you have 
used PyClaw for research, chances are that others would find your
code useful.  See :ref:`develop` for more details.


Citing
=======================
If you use PyClaw in work that is published, please cite the software::

    @misc{pyclaw,
    title={PyClaw software}, 
    url={http://numerics.kaust.edu.sa/pyclaw}, 
    author={Mandli, Kyle T and Ketcheson, David I. and others}, 
    note={Version x.y}
    year={2011}}

Please fill in the version number that you used.
As appropriate, please also cite any of the 
`publications based on PyClaw <http://www.mendeley.com/groups/1526933/pyclaw-publications/>`_.

If your work relies particularly on parallel capability, please cite::

    @inbook{petclaw2011,
    title={PetClaw: A Scalable Parallel Nonlinear Wave Propagation Solver for Python},
    booktitle={Proceedings of SpringSim 11}, 
    publisher={ACM},
    author={Alghamdi, Amal and Ahmadia, Aron and Ketcheson, David I. and Knepley, Matthew G. and Mandli, Kyle T and Dalcin, Lisandro}, 
    year={2011}}


License
=======================
PyClaw is distributed under a Berkeley Software Distribution
(BSD) style license.  The license is in the file pyclaw/LICENSE.txt and
reprinted below.

See http://www.opensource.org/licenses/bsd-license.php for more details.

Copyright (c) 2008-2011 Kyle Mandli and David I. Ketcheson.  All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright 
    notice, this list of conditions and the following disclaimer in the 
    documentation and/or other materials provided with the distribution.
  * Neither the name of King Abdullah University of Science & Technology nor 
    the names of its contributors may be used to endorse or promote products 
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Funding
==========

PyClaw development has been supported by 
grants from King Abdullah University of Science & Technology.
