#!/usr/bin/env python

# How to use this file
# python setup.py build_ext -i

import os
from os.path import join as pjoin
from os.path import pardir as pardir

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('euler_1d', parent_package, top_path)

    if top_path:
        # installing into clawpack
        sharpclaw_dir = pjoin(top_path, 'pyclaw', 'src', 'pyclaw', 'sharpclaw')
    else:
        # building as part of repository checkout
        this_dir = os.path.dirname(os.path.realpath(__file__))
        sharpclaw_dir = pjoin(this_dir, pardir, pardir, 'src', 'pyclaw', 'sharpclaw')
    sharpclaw_srcs = [pjoin(sharpclaw_dir, src) for src in ['ClawParams.f90',
                                                            'workspace.f90',
                                                            'weno.f90',
                                                            'reconstruct.f90','flux1.f90']]

    config.add_extension('sharpclaw1',
                         ['evec.f90'] + sharpclaw_srcs)
    config.add_extension('euler_tfluct','euler_tfluct.f90')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
