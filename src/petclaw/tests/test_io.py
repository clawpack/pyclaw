from __future__ import absolute_import
from clawpack import pyclaw
from clawpack import petclaw
import os

class PetClawIOTest(pyclaw.IOTest):
    @property
    def solution(self):
        return petclaw.Solution()

    @property
    def file_formats(self):
        return ['hdf5']

    @property
    def this_dir(self):
        return os.path.dirname(os.path.abspath(__file__))

    @property
    def test_data_dir(self):
        return os.path.join(self.this_dir, '../../pyclaw/tests/test_data')

    def test_io_from_binary(self):
        return