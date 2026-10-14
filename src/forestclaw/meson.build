package = 'clawpack.forestclaw'
pkg_dir = join_paths(package.split('.'))

python_sources = {
  '': [
  '__init__.py',
  'geometry.py',
  'test.py'
  ],
  'fileio': [
    '__init__.py',
    'ascii.py',
    'test.py'
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

