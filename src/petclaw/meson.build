package = 'clawpack.petclaw'
pkg_dir = join_paths(package.split('.'))

python_sources = {
  '': [
    '__init__.py',
    'cfl.py',
    'controller.py',
    'geometry.py',
    'plot.py',
    'solution.py',
    'state.py',
  ],
  'classic': [
    '__init__.py',
    'solver.py',
  ],
  'sharpclaw': [
    '__init__.py',
    'solver.py',
  ],
  'limiters':[
    '__init__.py',
  ],
  'fileio': [
    '__init__.py',
    'hdf5.py',
    'petsc.py',
  ]
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
