#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import sys

try:
    import numpy
except:
    sys.path.append("/opt/share/ksl/numpy/dev-aug29/ppc450d/lib/python/")
    import numpy

import ksl_rescale 

a=numpy.array([[1,2],[3,4]],dtype=float,order='FORTRAN')
print(a)
print("rescaling by 2!")
ksl_rescale.rescale(a,2.0)
print(a)
