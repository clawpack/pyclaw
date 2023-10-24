package = 'clawpack.pyclaw'
pkg_dir = join_paths(package.split('.'))

python_sources = {
  '': [
    '__init__.py',
    'cfl.py',
    'controller.py',
    'gauges.py',
    'geometry.py',
    'plot.py',
    'solution.py',
    'solver.py',
    'state.py',
    'util.py',
  ],
  'classic': [
    '__init__.py',
    'solver.py',
  ],
  'fileio': [
    '__init__.py',
    'ascii.py',
    'binary.py',
    'claw_vtk.py',
    'hdf5.py',
    'netcdf.py',
  ],
  'limiters': [
    '__init__.py',
    'recon.py',
    'reconstruct.py',
    'tvd.py',
  ],
  'limiters/weno': [
    '__init__.py',
  ],
  'sharpclaw': [
    '__init__.py',
    'solver.py',
    'weno.py',
  ],
  'tests': [
    '__init__.py',
    'test_io.py',
  ],
}

foreach subpkg, subsrcs: python_sources
  sources = []
  foreach src: subsrcs
    sources += subpkg / src
  endforeach
  py.install_sources(
    sources,
    subdir: pkg_dir / subpkg,
  )
endforeach


fortran_extensions = {
  'classic': {
    'classic1': ['limiter.f90','philim.f90','step1.f90'],
    'classic2': ['limiter.f90','philim.f90','flux2.f90','step2ds.f90','step2.f90'],
    'classic3': ['limiter.f90','philim.f90','flux3.f90','step3ds.f90','step3.f90'],
  },
  'sharpclaw': {
    'sharpclaw1': ['ClawParams.f90','weno.f90','reconstruct.f90','evec.f90','workspace.f90','flux1.f90'],
    'sharpclaw2': ['ClawParams.f90','weno.f90','reconstruct.f90','evec.f90','workspace.f90','flux2.f90','flux1.f90'],
    'sharpclaw3': ['ClawParams.f90','weno.f90','reconstruct.f90','evec.f90','workspace.f90','flux3.f90','flux1.f90'],
  },
}

foreach subpkg, module: fortran_extensions
  ext_srcs = []
  srcdir = join_paths(subpkg.split('.'))
  foreach ext_name, sources : module
    foreach src: sources
      ext_srcs += srcdir / src
    endforeach
    py.extension_module(
      ext_name, ext_srcs,
      incdir_f2py / 'fortranobject.c',
      include_directories: incdir_numpy,
      dependencies : py_dep,
      subdir: 'clawpack/' + srcdir,
      install : true
    )
  endforeach
endforeach
