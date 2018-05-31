[![Build Status](https://travis-ci.org/clawpack/pyclaw.svg?branch=master)](https://travis-ci.org/clawpack/pyclaw)
[![Coverage Status](https://img.shields.io/coveralls/clawpack/pyclaw.svg)](https://coveralls.io/r/clawpack/pyclaw?branch=master)

[![PyPI version](https://badge.fury.io/py/clawpack.svg)](https://badge.fury.io/py/clawpack)


Quick start:

    git clone git@github.com:clawpack/clawpack.git
    cd clawpack
    pip install -e .
    cd clawpack/pyclaw/examples/euler_2d
    python shockbubble.py iplot=1


# PyClaw

Pyclaw is a Python-based solver for hyperbolic PDEs that includes the algorithms
of Clawpack and SharpClaw.  
It has been designed with easy extensibility, performance, and exploration in mind.
PyClaw also includes a scalable parallel implementation of Clawpack using PETSc.

Full documentation is available at http://clawpack.org/pyclaw/.

You can get the latest development version of PyClaw from
http://github.com/clawpack/.

If you use PyClaw in research that leads to publication, please
[cite PyClaw](http://www.clawpack.org/pyclaw/index.html#citing).
