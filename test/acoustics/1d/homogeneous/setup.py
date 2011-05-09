#!/usr/bin/env python

#$ python setup.py build_ext --inplace

from numpy.distutils.command import build_src
from distutils.errors import DistutilsError

# a bit of monkeypatching ...
#import Cython.Compiler.Main
#build_src.build_src.Pyrex = Cython
#build_src.have_pyrex = lambda: True

def generate_a_cython_source(self, base, ext_name, source, extension):
    import os
    from distutils.dep_util import newer_group
    from numpy.distutils import log
    if self.inplace or not have_pyrex():
        target_dir = os.path.dirname(base)
    else:
        target_dir = appendpath(self.build_src, os.path.dirname(base))
    target_file = os.path.join(target_dir, ext_name + '.c')
    depends = [source] + extension.depends
    if self.force or newer_group(depends, target_file, 'newer'):
        import Cython.Compiler.Main
        log.info("pyrexc:> %s" % (target_file))
        self.mkpath(target_dir)
        options = Cython.Compiler.Main.CompilationOptions(defaults=Cython.Compiler.Main.default_options, include_path=extension.include_dirs, output_file=target_file)
        cython_result = Cython.Compiler.Main.compile(source, options=options)
        if cython_result.num_errors != 0:
            raise DistutilsError("%d errors while compiling %r with Cython" % (cython_result.num_errors, source))
        #elif os.path.isfile(target_file):
        #    log.warn("Cython required for compiling %r but not available,"\
        #                 " using old target %r"\
        #                 % (source, target_file))
    return target_file
build_src.build_src.generate_a_pyrex_source = generate_a_cython_source

def configuration(parent_package='',top_path=None):
    INCLUDE_DIRS = []
    LIBRARY_DIRS = []
    LIBRARIES    = []

    # PETSc
    import os
    PETSC_DIR  = os.environ['PETSC_DIR']
    PETSC_ARCH = os.environ.get('PETSC_ARCH', '')
    from os.path import join, isdir
    if PETSC_ARCH and isdir(join(PETSC_DIR, PETSC_ARCH)):
        INCLUDE_DIRS += [join(PETSC_DIR, PETSC_ARCH, 'include'),
                         join(PETSC_DIR, 'include')]
        LIBRARY_DIRS += [join(PETSC_DIR, PETSC_ARCH, 'lib')]
    else:
        if PETSC_ARCH: pass # XXX should warn ...
        INCLUDE_DIRS += [join(PETSC_DIR, 'include')]
        LIBRARY_DIRS += [join(PETSC_DIR, 'lib')]
    def haslibrary(basename):
        return (   os.path.exists(join(PETSC_DIR,PETSC_ARCH,'lib','%s.so' % basename))
                or os.path.exists(join(PETSC_DIR,PETSC_ARCH,'lib','%s.dylib' % basename)))
    if haslibrary('libpetsc'):
        LIBRARIES += ['petsc']
    elif haslibrary('libpetscts'):
        LIBRARIES += ['petscts','petscsnes','petscksp','petscdm','petscmat','petscvec','petscsys']
    else: raise DistutilsError('Cannot find valid PETSc libraries in PETSC_DIR=%s PETSC_ARCH=%s',PETSC_DIR,PETSC_ARCH)

    # PETSc for Python
    import petsc4py
    INCLUDE_DIRS += [petsc4py.get_include()]

    # Configuration
    from numpy.distutils.misc_util import Configuration
    module = 'assembly'
    config = Configuration('', parent_package, top_path)
    config.add_extension(module,
                         sources = [module+'.pyx',
                                    module+'impl.c'],
                         depends = [module+'impl.h'],
                         include_dirs=INCLUDE_DIRS + [os.curdir],
                         libraries=LIBRARIES,
                         library_dirs=LIBRARY_DIRS,
                         runtime_library_dirs=LIBRARY_DIRS)
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
