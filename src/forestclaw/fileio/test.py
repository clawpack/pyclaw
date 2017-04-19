#!/usr/bin/env python
# encoding: utf-8
r"""
Test suite for forestclaw.io
"""

from __future__ import absolute_import
from __future__ import print_function

import os
import tempfile
import shutil

import numpy
import nose

import clawpack.forestclaw as forestclaw


def test_forestclaw_input():
    """Simple test to read in a ForestClaw ASCII file"""

    # Create test solution
    x = forestclaw.geometry.Dimension(0.0, 1.0, 100, name='x')
    domain = forestclaw.geometry.Domain(x)
    state = forestclaw.state.State(domain, 1)
    state.q = numpy.zeros(x.num_cells)
    sol = forestclaw.solution.Solution(state, domain)

    # Test specific extensions to Patch
    sol.domain.patches[0].block_number = 2
    sol.domain.patches[0].mpi_rank = 2

    # Create temporary directory to write to and read from
    try:
        temp_path = tempfile.mkdtemp()

        # Real test
        sol.write(0, path=temp_path)
        read_sol = forestclaw.Solution(0, path=temp_path)
        nose.tools.assert_equal(sol, read_sol,
                                "ForestClaw IO read/write test failed, " +
                                "solutions are not equal.")
    # except:
        # shutil.copytree(temp_path, os.getcwd())
    finally:
        shutil.rmtree(temp_path)

if __name__ == '__main__':
    import nose
    nose.main()
