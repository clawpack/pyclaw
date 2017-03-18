#!/usr/bin/env python
# encoding: utf-8
r"""
Test suite for forestclaw.io
"""

import os
import tempfile

import nose

import clawpack.forestclaw as forestclaw

def test_forestclaw_input():
    """Simple test to read in a ForestClaw ASCII file"""

    # Create test solution
    x = forestclaw.geometry.Dimension(0.0, 1.0, 100, name='x')
    state = forestclaw.state(1, 1)
    state.q = numpy.zeros(x.num_cells)
    sol = forestclaw.Solution(forestclaw.geometry.Domain(x), state)

    # Test specific extensions to Patch
    sol.patches[0].block_number = 2
    sol.patches[0].mpi_rank = 2

    # Create temporary directory to write to and read from
    try:
        temp_path = tempfile.mkdtemp()

        # Real test
        sol.write(0, path=temp_path)
        read_sol = forestclaw.Solution(0, path=temp_path)
        nose.tools.assert_equal(sol, read_sol,
                                "ForestClaw IO read/write test failed, " +
                                "solutions are not equal.")
    except:
        shutil.copytree(self.temp_path, os.getcwd())
    shutil.rmtree(self.temp_path)

if __name__ == '__main__':
    import nose
    nose.main()
