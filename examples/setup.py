#!/usr/bin/env python
import os
import os.path

def configuration(parent_package='',top_path=None):
    """Configuration for examples directory

    Since we don't want to manually maintain the examples setup file,
    this setup file treats all subdirectories as Python packages.  The
    developer is still responsible for putting an __init__.py file in
    each subdirectory.
    """

    from numpy.distutils.misc_util import Configuration
    config = Configuration('examples', parent_package, top_path)

    # automatic detection of subdirectories in path
    d = os.path.dirname(os.path.realpath(__file__))
    subdirs = (s for s in os.listdir(d) if os.path.isdir(os.path.join(d,s)))
    for package in subdirs:
        config.add_subpackage(package)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
