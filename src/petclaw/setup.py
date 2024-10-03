#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('petclaw', parent_package, top_path)
    config.add_data_files('log.config')
    config.add_subpackage('classic')
    config.add_subpackage('sharpclaw')
    config.add_subpackage('limiters')
    config.add_subpackage('fileio')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
