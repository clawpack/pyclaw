package = 'clawpack.pyclaw'
pkg_dir = join_paths(package.split('.'))

python_sources = [
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
  'classic/__init__.py',
  'classic/solver.py',
  'sharpclaw/__init__.py',
  'sharpclaw/solver.py',
  'sharpclaw/weno.py',
  'fileio/__init__.py',
  'fileio/ascii.py',
  'fileio/binary.py',
  'fileio/claw_vtk.py',
  'fileio/hdf5.py',
  'fileio/netcdf.py',
  'limiters/__init__.py',
  'limiters/recon.py',
  'limiters/reconstruct.py',
  'limiters/tvd.py',
  'limiters/weno/__init__.py',
  'tests/__init__.py',
  'tests/test_io.py',
]

py.install_sources(
  python_sources,
  subdir: pkg_dir,
)

