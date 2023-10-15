# PyClaw examples

These directories contain several prepared examples of how to run PyClaw.
There are multiple ways to run and plot an example:

### Method 1: From the command line

    cd examples/acoustics_1d_homogeneous
    python acoustics_1d.py iplot=1

The command-line option `iplot=1` will cause an interactive plotting session
to be launched after the simulation is complete.

### Method 2: From IPython

    ipython --matplotlib
    from clawpack.pyclaw import examples
    claw = examples.acoustics_1d.setup()
    claw.run()
    claw.plot()

Some of the more basic examples are the homogeneous acoustics, advection,
burgers, traffic, and shallow water examples.  Advanced examples include
`Rossby_wave` (2D grid mapped to a sphere) `shock_bubble_interaction` (cylindrical
symmetry with geometric source terms), and `psystem_2d` (shows how to set up gauges, 
a `before_step` function, and more).
