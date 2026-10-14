from clawpack import pyclaw
from clawpack.pyclaw.solution import Solution

class Solution(Solution):
    """ Parallel Solution class.
    """
    __doc__ += pyclaw.util.add_parent_doc(pyclaw.Solution)

    def get_read_func(self, file_format):
        from clawpack.petclaw import fileio

        if file_format == 'petsc':
            try:
                import clawpack.petclaw.fileio
                return clawpack.petclaw.fileio.petsc.read
            except AttributeError as e:
                try:
                    from petsc4py import PETSc
                except ImportError:
                    raise ImportError("petsc4py is required for reading petsc format files, but petsc4py could not be imported.")
                raise e
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
