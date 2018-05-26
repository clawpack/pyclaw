#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing an ascii output file from ForestClaw
"""

from __future__ import absolute_import
from __future__ import print_function

import logging
import numpy as np
# import clawpack.forestclaw as forestclaw

from .. import Dimension, Patch
from clawpack.pyclaw.fileio.ascii import write
from clawpack.pyclaw.fileio.ascii import read
from clawpack.pyclaw.util import read_data_line
import forestclaw



logger = logging.getLogger('pyclaw.fileio')


def write_patch_header(f, patch):
    f.write("%5i                  patch_number\n" % patch.patch_index)
    f.write("%5i                  AMR_level\n" % patch.level)
    f.write("%5i                  block_number\n" % patch.block_number)
    f.write("%5i                  mpi_rank\n" % patch.mpi_rank)
    for dim in patch.dimensions:
        f.write("%5i                  m%s\n" % (dim.num_cells, dim.name))
    for dim in patch.dimensions:
        f.write("%18.8e     %slow\n" % (dim.lower, dim.name))
    for dim in patch.dimensions:
        f.write("%18.8e     d%s\n" % (dim.delta, dim.name))

    f.write("\n")


def read_patch_header(f, num_dim):
    r"""Read header describing the next patch

    :Input:
     - *f* - (file) Handle to open file
     - *num_dim* - (int) Number of dimensions

    :Output:
     - *patch* - (clawpack.pyclaw.geometry.Patch) Initialized patch represented
       by the header data.
    """

    n = np.zeros((num_dim))
    d = np.zeros((num_dim))
    lower = np.zeros((num_dim))

    patch_index = read_data_line(f, data_type=int)
    level = read_data_line(f, data_type=int)
    block_number = read_data_line(f, data_type=int)
    mpi_rank = read_data_line(f, data_type=int)
    for i in range(num_dim):
        n[i] = read_data_line(f, data_type=int)
    for i in range(num_dim):
        lower[i] = read_data_line(f)
    for i in range(num_dim):
        d[i] = read_data_line(f)

    blank = f.readline()

    # Construct the patch
    # Since we do not have names here, we will construct the patch with
    # dimension names x,y,z
    names = ['x', 'y', 'z']
    dimensions = [forestclaw.geometry.Dimension(lower[i], lower[i] + n[i] * d[i],
                  n[i], name=names[i]) for i in range(num_dim)]
    patch = forestclaw.geometry.Patch(dimensions)

    # Add AMR attributes:
    patch.patch_index = patch_index
    patch.level = level

    # Add ForestClaw attributes
    patch.block_number = block_number
    patch.mpi_rank = mpi_rank

    return patch

# # Replace the ascii module functions with those defined above
# # ascii.read_patch_header = read_patch_header
# # ascii.write_patch_header = write_patch_header
# read = ascii.read
# write = ascii.write
