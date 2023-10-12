package = 'clawpack.petclaw'
pkg_dir = join_paths(package.split('.'))

python_sources = [
  '__init__.py',
  'cfl.py',
  'controller.py',
  'geometry.py',
  'plot.py',
  'solution.py',
  'state.py',
  'fileio/__init__.py',
  'fileio/hdf5.py',
  'fileio/petsc.py',
  'limiters/__init__.py',
  'tests/test_io.py',
  'classic/__init__.py',
  'classic/solver.py',
  'sharpclaw/__init__.py',
  'sharpclaw/solver.py',
  'tests/test_io.py',
]

py.install_sources(
  python_sources,
  subdir: pkg_dir,
)

