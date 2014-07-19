#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('sharpclaw', parent_package, top_path)

    config.add_extension('sharpclaw1',
                         ['ClawParams.f90','weno.f90','reconstruct.f90',
                          'evec.f90','workspace.f90','flux1.f90'])

    config.add_extension('sharpclaw2',
                         ['ClawParams.f90','weno.f90','reconstruct.f90',
                          'evec.f90','workspace.f90','flux2.f90',
                          'flux1.f90'])

    config.add_extension('sharpclaw3',
                         ['ClawParams.f90','weno.f90','reconstruct.f90',
                          'evec.f90','workspace.f90','flux3.f90',
                          'flux1.f90'])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
