This folder contains the following files

Makefile: to build the FRS.so module from the fortran code main.90 created by Kyle Mandli in ../ directory. You need to have f2py tool and gfortran from http://r.research.att.com/tools/. Hpc gfortran 4.5 version did not work for me. 

RiemannSolver.py: is the python wrapper for main.90. It imports FRS and contains one class that provides several methods.

riemannSolverDemo.py: is a demonstration example of using class RiemannSolver.py and plotting some of the results

riemannSolverAssertion.py: is a test file to assert the result of computation with a certain tolerance. It reads the input from q.npy and aux.npy, calculate the output and compare it with a previously calculated output stored in s.npy and waves.npy. (it needs to be modified to read all the inputs not only q and aux)
 
riemannSolverSerialTimer.py: This script calculates the execution time for different input sizes and different numbers of time steps and plots the timing results.
