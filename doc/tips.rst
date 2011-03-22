Some tips
=======================

PetClaw now uses a different array indexing than Clawpack.  In PetClaw,
the value of the $m$th conserved quantity at $(x_i,y_j)$ is ::

    q[m,i,j]

That is, the index $m$ comes first, whereas in Clawpack it comes last.
This "interleaved" array layout is potentially more cache-efficient.
The next version of Clawpack (5.0) will also use interleaved arrays.

When using the PyClaw controller, it is necessary to set ::

        controller.output_format = 'petsc'

for compatibility with PetClaw.  Note that this is not the default.
Similarly, when reading in data for analysis or plotting, it is necessary
to set the format to 'petsc'.
