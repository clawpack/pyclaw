from clawpack import pyclaw
from clawpack import petclaw
import os

class TestParallelIO(pyclaw.TestIO):
    @property
    def solution(self):
        return petclaw.Solution()

    @property
    def file_formats(self):
        return ['hdf5']

    def test_io_from_binary(self):
        return
