.. _clawpack_and_pyclaw:

======================================================
Clawpack and PyClaw
======================================================

.. _port_Example:

Porting a problem from Clawpack to PyClaw 
======================================================
The script `pyclaw/development/clawdata2pyclaw` is intended to aid
in converting a Clawpack 4.6 problem setup to PyClaw.  However,
some manual conversion is necessary, especially if the problem
includes custom fortran routines.

In PyClaw, the high-level portions of the Fortran routines are reorganized in 
an object-oriented Python framework, while the low-level ones are bound through
the Fortran to Python interface generator `f2py <http://www.scipy.org/F2py>`_.
Therefore, in a typical situation the user is shielded from f2py. However, if 
the user wants to reutilize some problem-specific fortran routines that were set up and 
tested in a Clawpack problem, he can easily do it. Indeed, if those routines 
are complicated and implement time consuming algorithms, and the performance 
of the overall procedure (porting and running into PyClaw a new problem) is one 
of the main user's concern, one should consider directly using the f2py 
interface in the initialization script (see :ref:`problem_setup`).
The shallow water equations solved on a sphere `(code here) <http://numerics.kaust.edu.sa/pyclaw/apps/shallow-sphere/shallow_4_Rossby_Haurwitz_wave.py>`_ represent a
useful and complete example to understand the simplicity of the procedure. 
In that problem setup, a few Fortran routines have been used to provide the 
following functionality:

    * Initialize the solution ``state.q[:,:,:]``

    * Provide the mapping from a uniform Cartesian domain to the desired 
      physical domain, i.e. the mapc2p function

    * Setup the auxiliary variables ``state.aux[:,:,:]``

    * Compute the (non-hyperbolic) contribution of a source term

    * Impose custom boundary conditions to both solution and auxiliary 
      variables

The first step to succesfully interface the Fortran functions with PyClaw 
is to automate the extension module generation of these routines through f2py.
`This Makefile
<http://numerics.kaust.edu.sa/pyclaw/apps/shallow-sphere/shallow_4_Rossby_Haurwitz_wave.py>`_
shows how to do it::

    # Problem's source Fortran files
    INITIALIZE_SOURCE = mapc2p.f setaux.f qinit.f src2.f
    problem.so: $(INITIALIZE_SOURCE)
        $(F2PY) -m problem -c $^

In the code above, we are giving to f2py the instructions to compile a 
set of Fortran routines (the INITIALIZE_SOURCE container) and build a module 
(``problem.so``) which can then be imported into Python and used there like a normal
function. Indeed, f2py scans Fortran codes to produce the signature files (.pyf files)
which contain all the information (function names, arguments and 
their types, etc.) that is needed to construct Python bindings to Fortran 
functions. The argument following the ''-m'' flag is the name the python module should have (i.e.
the name of the target). f2py uses the ``numpy.distutils`` module from NumPy 
that supports a number of major Fortran compilers. For more information please 
look at `<http://www.scipy.org/F2py>`_.

After the compilation has been succesfully completed, the signature of each 
function contained in ``problem.so`` must be verified and the intent of the 
variables added (if there was nothing stated in the 
code). One can easily achieve that by using the following commands::

    $ ipython
    >>> import problem
    >>> problem?

The last command queries the content of the module and outputs the functions' 
signature that must be used in the initialization script to correctly call the 
fortran functions. In the shallow water equations on a sphere example, we get 
the following output::

    >>> Type:		module
    >>> Base Class:	<type 'module'>
    >>> String Form:	<module 'problem' from 'problem.so'>
    >>> Namespace:	Interactive
    >>> File:		/Users/../../../clawpack/pyclaw/apps/shallow-sphere/problem.so
    >>> Docstring:
        This module 'problem' is auto-generated with f2py (version:1).
        Functions:
        mapc2p(x1,y1,xp,yp,zp,rsphere)
        aux = setaux(maxmx,maxmy,num_ghost,mx,my,xlower,ylower,dxc,dyc,aux,rsphere,num_aux=shape(aux,0))
        q = qinit(maxmx,maxmy,num_ghost,mx,my,xlower,ylower,dx,dy,q,aux,rsphere,num_eqn=shape(q,0),num_aux=shape(aux,0))
        q = src2(maxmx,maxmy,num_ghost,xlower,ylower,dx,dy,q,aux,t,dt,rsphere,num_eqn=shape(q,0),mx=shape(q,1),my=shape(q,2),num_aux=shape(aux,0))

For instance, the function ``src2``, which computes the contribution of the 
(non-hyperbolic) source term, has the following intent variables::

    >>> cf2py integer intent(in) maxmx
    >>> cf2py integer intent(in) maxmy
    >>> cf2py integer optional, intent(in) num_eqn
    >>> cf2py integer intent(in) num_ghost
    >>> cf2py integer intent(in) mx
    >>> cf2py integer intent(in) my
    >>> cf2py double precision intent(in) xlower
    >>> cf2py double precision intent(in) ylower
    >>> cf2py double precision intent(in) dx
    >>> cf2py double precision intent(in) dy
    >>> cf2py intent(in,out) q
    >>> cf2py integer optional, intent(in)  num_aux
    >>> cf2py intent(in) aux
    >>> cf2py double precision intent(in) t
    >>> cf2py double precision intent(in) dt
    >>> cf2py double precision intent(in) Rsphere

Note that ``num_eqn``, ``mx``, ``my`` ``num_aux`` are identified by f2py as optional
arguments since their values can be retrieved by looking at the dimensions of
other multidimensional arrays, i.e. ``q`` and ``aux``.

We are now ready to call and use the Fortran functions in the initialization
script. For instance, the ``src2`` function is called in the 
`script <http://numerics.kaust.edu.sa/pyclaw/apps/shallow-sphere/shallow_4_Rossby_Haurwitz_wave.py>`_ by using a fortran_src_wrapper function whose main part reads::

    >>> # Call src2 function
    >>> import problem
    >>> state.q = problem.src2(mx,my,num_ghost,xlowerg,ylowerg,dx,dy,q,aux,t,dt,Rsphere)

A similar approach is used to call other wrapped Fortran functions like 
``qinit``, ``setaux``, etc.

An important feature that makes PyClaw more flexible is the 
capability to replace the standard low-level Fortran routines whith some 
problem-specific routines. Binding new low-level functions and replacing the 
standard ones is very easy; the user just needs to modify the problem-specific 
Makefile. The shallow water equations on a sphere is again a 
typical example that uses this nice feature. Indeed, to run correctly the problem an 
ad-hoc ``step2`` function (i.e. the ``step2qcor``) is required. For that problem
the interesting part of the `Makefile
<http://numerics.kaust.edu.sa/pyclaw/apps/shallow-sphere/shallow_4_Rossby_Haurwitz_wave.py>`_
reads::

    # Override step2.f with a new function that contains a call to an additional
    # function, i.e. qcor.f
    # ==========================================================================
    override TWO_D_CLASSIC_SOURCES = step2qcor.f qcor.o flux2.o limiter.o philim.o

    qcor.o: qcor.f
        $(FC) $(FFLAGS) -o qcor.o -c qcor.f

The user has just to override ``step2.f`` with the new function ``step2qcor.f`` 
and provide new::

    output_filenames : input_filenames
    	actions

rules to create the targets required by the new Fortran routine. 
Similar changes to the problem-specific Makefile can be used to replace other 
low-level Fortran routines.


.. _diffs:

Important differences between PyClaw and Clawpack
==================================================
PyClaw incorporates some important changes relative to Clawpack.  
Most of these are planned for inclusion in Clawpack 5.0.

Interleaved arrays
===================
PyClaw uses a different array indexing than Clawpack.  In PetClaw,
the value of the :math:`m`-th conserved quantity at :math:`(x_i,y_j)` is ::

    q[m,i,j]

That is, the index :math:`m` comes first, whereas in Clawpack it comes last.
This "interleaved" array layout is potentially more cache-efficient.
The next version of Clawpack (5.0) will also use interleaved arrays.

