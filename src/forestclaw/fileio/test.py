#!/usr/bin/env python
# encoding: utf-8
r"""
Test suite for forestclaw.io
"""

import os
import shutil
import tempfile
import shutil

import numpy
import nose

from .. import geometry
from clawpack.pyclaw.solution import Solution
from clawpack.pyclaw.state import State


def test_forestclaw_input():
    """Simple test to read in a ForestClaw ASCII file"""

    # Create test solution
    x = geometry.Dimension(0.0, 1.0, 100, name='x')
    state = State(geometry.Patch(x), 1)
    state.q = numpy.zeros((1, x.num_cells))
    sol = Solution(state, geometry.Domain(x))

    # Test specific extensions to Patch
    sol.domain.patches[0].block_number = 2
    sol.domain.patches[0].mpi_rank = 2

    # Create temporary directory to write to and read from
    try:
        temp_path = tempfile.mkdtemp()

        # Real test
        sol.write(0, path=temp_path)
        read_sol = Solution(0, path=temp_path,
                                       file_format='forestclaw')
        nose.tools.assert_equal(sol, read_sol,
                                "ForestClaw IO read/write test failed, " +
                                "solutions are not equal.")
    except:
        shutil.copy(os.path.join(temp_path, 'fort.q0000'), os.getcwd())
    finally:
        shutil.rmtree(temp_path)

if __name__ == '__main__':
    nose.main()
