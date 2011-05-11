Some tips
=======================

When using the PyClaw controller, it is necessary to set ::

        controller.output_format = 'petsc'

for compatibility with PetClaw.  Note that this is not the default.
Similarly, when reading in data for analysis or plotting, it is necessary
to set the format to 'petsc'.
