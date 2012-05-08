r"""Distutils script to build the WENO reconstructor."""

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('weno', parent_package, top_path)

    config.add_extension('reconstruct',
                         sources = ['reconstruct.c'],
                         include_dirs = [get_numpy_include_dirs()],
                         extra_compile_args = ['-std=c99'])    # to support 'restrict' with gcc                               
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
