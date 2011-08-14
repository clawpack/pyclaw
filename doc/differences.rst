.. _diffs:

==================================================
Important differences between PyClaw and Clawpack
==================================================
PyClaw incorporates some important changes relative to Clawpack.  
Most of these are planned for inclusion in Clawpack 5.0.

Interleaved arrays
===================
PyClaw uses a different array indexing than Clawpack.  In PetClaw,
the value of the :math:`m`-th conserved quantity at :math:`(x_i,y_j)` is ::

    q[m,i,j]

That is, the index $m$ comes first, whereas in Clawpack it comes last.
This "interleaved" array layout is potentially more cache-efficient.
The next version of Clawpack (5.0) will also use interleaved arrays.


