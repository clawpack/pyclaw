#!/usr/bin/env python

import os

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('classic', parent_package, top_path)

    config.add_extension('classic1',
                         ['limiter.f90','philim.f90','step1.f90'],f2py_options=['--quiet'])

    config.add_extension('classic2',
                         ['limiter.f90','philim.f90','flux2.f90','step2ds.f90','step2.f90'],f2py_options=['--quiet'])

    config.add_extension('classic3',
                         ['limiter.f90','philim.f90','flux3.f90','step3ds.f90','step3.f90'],f2py_options=['--quiet'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
