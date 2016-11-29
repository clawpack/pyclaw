#!/usr/bin/env python
# encoding: utf-8
r"""
This script is for testing petsc4py setup. If you
run this script with four processes as follows: ::

    $ mpiexec -n 4 python petsc_hello_world.py 

Then the expected output should look like the following: ::
        
    Hello World! From process 3 out of 4 process(es).
    Hello World! From process 1 out of 4 process(es).
    Hello World! From process 0 out of 4 process(es).
    Hello World! From process 2 out of 4 process(es).

"""
from __future__ import absolute_import
from __future__ import print_function
from petsc4py import PETSc

rank = PETSc.COMM_WORLD.getRank()
size = PETSc.COMM_WORLD.getSize()

print('Hello World! From process {rank} out of {size} process(es).'.format(rank=rank,size=size))
