from __future__ import absolute_import
from clawpack import pyclaw
from clawpack.pyclaw.solution import Solution

class Solution(Solution):
    """ Parallel Solution class.
    """
    __doc__ += pyclaw.util.add_parent_doc(pyclaw.Solution)

    def get_read_func(self, file_format):
        from clawpack.petclaw import fileio
        if file_format == 'petsc':
            return fileio.petsc.read
        elif file_format == 'hdf5':
            return fileio.hdf5.read
        else:
            raise ValueError("File format %s not supported." % file_format)

    def get_write_func(self, file_format):
        from clawpack.petclaw import fileio
        if 'petsc' in file_format:
            return fileio.petsc.write
        elif 'hdf5' in file_format:
            return fileio.hdf5.write
        else:
            raise ValueError("File format %s not supported." % file_format)


