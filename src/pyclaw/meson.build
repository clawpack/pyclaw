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

examples = {
  '' : [
    '__init__.py'
  ],
  'acoustics_1d_homogeneous' : [
    '__init__.py',
    'acoustics_1d.py',
    'test_acoustics.py',
  ],
  'acoustics_2d_homogeneous' : [
    '__init__.py',
    'acoustics_2d.py',
    'test_2d_acoustics.py',
  ],
  'acoustics_2d_mapped' : [
    '__init__.py',
    'acoustics_2d_inclusions.py',
    'test_acoustics_2d_mapped.py',
  ],
  'acoustics_2d_variable' : [
    '__init__.py',
    'acoustics_2d_interface.py',
    'test_acoustics_2d_variable.py',
    'test_acoustics_2d_variable_io.py',
  ],
  'acoustics_3d_variable' : [
    '__init__.py',
    'acoustics_3d_interface.py',
    'test_3d_acoustics.py',
  ],
  'advection_1d' : [
    '__init__.py',
    'advection_1d.py',
    'advection_1d_nonunif.py',
    'test_advection.py',
    'test_advection_nonunif.py',
  ],
  'advection_1d_variable' : [
    '__init__.py',
    'variable_coefficient_advection.py',
    'test_variable_coefficient_advection.py',
  ],
  'advection_2d' : [
    '__init__.py',
    'advection_2d.py',
    'test_2d_advection.py',
  ],
  'advection_2d_annulus' : [
    '__init__.py',
    'advection_annulus.py',
    'test_2d_advection_annulus.py',
    'mapc2p.py',
  ],
  'advection_reaction_2d' : [
    '__init__.py',
    'advection_reaction.py',
    'test_advection_reaction.py',
  ],
  'burgers_1d' : [
    '__init__.py',
    'burgers_1d.py',
    'test_burgers_1d.py',
  ],
  'cubic_1d' : [
    '__init__.py',
    'cubic.py',
    'test_cubic.py',
  ],
  'euler_1d' : [
    '__init__.py',
    'shocksine.py',
    'shocktube.py',
    'test_shocksine.py',
    'test_shocktube.py',
    'test_woodward_colella_blast.py',
    'woodward_colella_blast.py',
  ],
  'euler_2d' : [
    '__init__.py',
    'euler_2d.py',
    'quadrants.py',
    'shock_bubble_interaction.py',
    'shock_forward_step.py',
    'test_euler2d.py',
  ],
  'euler_3d' : [
    'Sedov.py',
    '__init__.py',
    'shock_bubble.py',
    'shocktube.py',
    'test_euler3d.py',
  ],
  'euler_gravity_3d' : [
    '__init__.py',
    'plotCreateVisitXDMF.py',
    'rising_hot_sphere.py',
    'rising_hot_sphere_spherical.py',
    'test_eulerg3d.py',
    'visitHDF5XDMF.py',
  ],
  'iso_c_advection' : [
    '__init__.py',
    'iso_c_advection.py',
    'iso_c_solver.py',
  ],
  'kpp' : [
    '__init__.py',
    'kpp.py',
  ],
  'mhd_1d' : [
    '__init__.py',
    'shocktube.py',
    'test_shocktube.py',
  ],
  'peano_shallow_2d' : [
    'setplot.py',
    'shallow2D.py',
    'test_identical_grids.py',
    'test_peano_solver.py',
  ],
  'psystem_2d' : [
    '__init__.py',
    'psystem_2d.py',
    'test_2d_psystem.py',
  ],
  'shallow_1d' : [
    '__init__.py',
    'dam_break.py',
    'sill.py',
    'test_shallow1d.py',
  ],
  'shallow_2d' : [
    '__init__.py',
    'radial_dam_break.py',
    'sill.py',
    'test_shallow2d.py',
  ],
  'shallow_sphere' : [
    '__init__.py',
    'setplot.py',
    'test_shallow_sphere.py',
    'Rossby_wave.py',
  ],
  'stegoton_1d' : [
    '__init__.py',
    'stegoton.py',
    'test_stego.py',
  ],
  'traffic' : [
    '__init__.py',
    'traffic.py',
    'test_traffic.py',
  ],
}

foreach ex, exsrcs: examples
  sources = []
  foreach src: exsrcs
    sources += '../../examples' / ex / src
  endforeach
  py.install_sources(
    sources,
    subdir: pkg_dir / 'examples' / ex
  )
endforeach

py.install_sources(
    'log.config',
    subdir: pkg_dir
)

fortran_extensions = {
  'classic': {
    'classic1': ['limiter.f90','philim.f90','step1.f90'],
    'classic2': ['limiter.f90','philim.f90','flux2.f90','step2ds.f90','step2.f90'],
    'classic3': ['limiter.f90','philim.f90','flux3.f90','step3ds.f90','step3.f90'],
    'classic2_sw_sphere': ['limiter.f90','philim.f90','flux2.f90','step2ds.f90','../../../examples/shallow_sphere/step2qcor.f90','../../../examples/shallow_sphere/qcor.f90'],
  },
  'sharpclaw': {
    'sharpclaw1': ['ClawParams.f90','weno.f90','reconstruct.f90','evec.f90','workspace.f90','flux1.f90'],
    'sharpclaw2': ['ClawParams.f90','weno.f90','reconstruct.f90','evec.f90','workspace.f90','flux2.f90','flux1.f90'],
    'sharpclaw3': ['ClawParams.f90','weno.f90','reconstruct.f90','evec.f90','workspace.f90','flux3.f90','flux1.f90'],
    # 1D Euler with characteristic decomposition for SharpClaw limiting
    'euler_sharpclaw1': ['ClawParams.f90','weno.f90','reconstruct.f90','../../../examples/euler_1d/evec.f90','workspace.f90','flux1.f90'],
    'euler_tfluct1': ['../../../examples/euler_1d/euler_tfluct.f90'],
  },
}

numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

foreach subpkg, module: fortran_extensions
  srcdir = join_paths(subpkg.split('.'))
  foreach ext_name, sources : module
    ext_srcs = []
    foreach src: sources
      ext_srcs += srcdir / src
    endforeach
    f2py_extra = []
    #if subpkg == 'sharpclaw'
    if 'ClawParams.f90' in sources
        f2py_extra += ext_name + '-f2pywrappers2.f90'
    endif
    f2py_srcs = custom_target(
      command: [f2py, ext_name],
      input: ext_srcs,
      output: [ext_name + 'module.c', ext_name + '-f2pywrappers.f', f2py_extra],
    )
    py.extension_module(
      ext_name, [ext_srcs, f2py_srcs],
      incdir_f2py / 'fortranobject.c',
      c_args: [numpy_nodepr_api],
      include_directories: inc_np,
      subdir: 'clawpack/pyclaw/' + srcdir,
      install : true
    )
  endforeach
endforeach

py.extension_module(
    'reconstruct', 
    ['limiters/weno/reconstruct.c'],
    include_directories: inc_np,
    subdir: 'clawpack/pyclaw/limiters/weno',
    install : true
)

# Special extensions for examples

## Shallow water on the sphere

ext_name = 'sw_sphere_problem'
ext_srcs = [
  '../../examples/shallow_sphere/mapc2p.f90',
  '../../examples/shallow_sphere/setaux.f90',
  '../../examples/shallow_sphere/qinit.f90',
  '../../examples/shallow_sphere/src2.f90',
]
f2py_srcs = custom_target(
  command: [f2py, ext_name],
  input: ext_srcs,
  output: [ext_name + 'module.c', ext_name + '-f2pywrappers.f'],
)

py.extension_module(
  ext_name, [ext_srcs, f2py_srcs],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  subdir: 'clawpack/pyclaw/examples/shallow_sphere',
  install : true
)

## 2D Advection-reaction problem

ext_name = 'advection_2d'
ext_srcs = [
  '../../examples/advection_reaction_2d/rpn2_vc_advection_vector.f90',
  '../../examples/advection_reaction_2d/rpt2_dummy.f90',
]
f2py_srcs = custom_target(
  command: [f2py, ext_name],
  input: ext_srcs,
  output: [ext_name + 'module.c', ext_name + '-f2pywrappers.f'],
)

py.extension_module(
  ext_name, [ext_srcs, f2py_srcs],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  subdir: 'clawpack/pyclaw/examples/advection_reaction_2d',
  install : true
)
