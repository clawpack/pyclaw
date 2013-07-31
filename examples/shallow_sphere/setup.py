#!/usr/bin/env python

# How to use this file
# python setup.py build_ext -i

import os

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)

    config.add_extension('classic2',
                         ['step2qcor.f90','qcor.f90','../../src/pyclaw/classic/limiter.f90','../../src/pyclaw/classic/philim.f90','../../src/pyclaw/classic/flux2.f90','../../src/pyclaw/classic/step2ds.f90'])

    config.add_extension('problem',
                         ['mapc2p.f90','setaux.f90','qinit.f90','src2.f90'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
