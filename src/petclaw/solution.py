from clawpack import pyclaw
from clawpack.pyclaw.solution import Solution

class Solution(Solution):
    """ Parallel Solution class.
    """
    __doc__ += pyclaw.util.add_parent_doc(pyclaw.Solution)
