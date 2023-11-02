def configuration(parent_package='',top_path=None):
    import os

    from numpy.distutils.misc_util import Configuration

    config = Configuration('iso_c_advection', parent_package, top_path)

    print("You will probably want to run make manually here")

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration().todict())
