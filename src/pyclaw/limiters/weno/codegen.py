r"""Generate C code to compute a 2k-1 order WENO reconstructions.

This generates a Python extension module (written in C) to perform
WENO reconstructions.

Usage:

$ python codegen.py
$ python setup.py build
$ cp build/lib*/reconstruct.so .

Note: the naming convection in PyWENO is cell based: 'left' means the
left edge of a cell, and 'right' means the right edge of a cell.
Matching this up with how PyClaw expects 'ql' and 'qr' to be indexed
is a bit tricky.

"""

from __future__ import absolute_import
from __future__ import print_function
import os

import pyweno.symbolic
import pyweno.c
from six.moves import range

# config

K = list(range(3, 10))                        # 2*k-1 order reconstructions
module = 'reconstruct'                  # py module name
output = module + '.c'

# open output and init

f = open(output, 'w')

c = pyweno.c.CCodeGenerator()
f.write(c.wrapper_head(module))

# smoothness
for k in K:

    print('generating code for k = %d...' % k)

    beta = pyweno.symbolic.jiang_shu_smoothness_coefficients(k)
    c.set_smoothness(beta)

    f.write(c.uniform_smoothness(function='smoothness_k'+str(k),
                                 wrapper=True))

    # left edge of cell (right side of boundary)

    (varpi, split) = pyweno.symbolic.optimal_weights(k, 'left')
    coeffs = pyweno.symbolic.reconstruction_coefficients(k, 'left')

    c.set_optimal_weights(varpi, split)
    c.set_reconstruction_coefficients(coeffs)

    f.write(c.uniform_weights(function='weights_left_k' + str(k),
                              wrapper=True))
    f.write(c.uniform_reconstruction(function='reconstruct_left_k' + str(k),
                                     wrapper=True))

    # right edge of cell (left side boundary, shifted)

    (varpi, split) = pyweno.symbolic.optimal_weights(k, 'right')
    coeffs = pyweno.symbolic.reconstruction_coefficients(k, 'right')

    c.set_optimal_weights(varpi, split)
    c.set_reconstruction_coefficients(coeffs)

    f.write(c.uniform_weights(function='weights_right_k' + str(k),
                              wrapper=True))
    f.write(c.uniform_reconstruction(function='reconstruct_right_k' + str(k),
                                     wrapper=True))


f.write(c.wrapper_foot())
f.close()

try:
    import os
    os.system('indent ' + output)
except:
    pass
