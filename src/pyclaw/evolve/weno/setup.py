r"""Distutils script to build the WENO reconstructor."""

import numpy as np
from distutils.core import setup, Extension

reconstruct = Extension('reconstruct',
                        sources = ['reconstruct.c'],
                        include_dirs = [np.get_include()],
                        extra_compile_args = ['-std=c99'],    # to support 'restrict' with gcc
                        )

setup(name = 'WENO reconstructor',
      version = '0.1',
      ext_modules = [reconstruct])
