#!/usr/bin/env python

# How to use this file
# python setup.py build_ext -i

from __future__ import absolute_import
import os
from os.path import join as pjoin
from os.path import pardir as pardir

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('shallow_sphere', parent_package, top_path)

    if top_path:
        # installing into clawpack
        classic_dir = pjoin(top_path, 'pyclaw', 'src', 'pyclaw', 'classic')
    else:
        # building as part of repository checkout
        this_dir = os.path.dirname(os.path.realpath(__file__))
        classic_dir = pjoin(this_dir, pardir, pardir, 'src', 'pyclaw', 'classic')
    classic_srcs = [pjoin(classic_dir, src) for src in ['limiter.f90',
                                                        'philim.f90',
                                                        'flux2.f90',
                                                        'step2ds.f90']]

    config.add_extension('classic2',
                         ['step2qcor.f90','qcor.f90'] + classic_srcs)

    config.add_extension('problem',
                         ['mapc2p.f90','setaux.f90','qinit.f90','src2.f90'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
