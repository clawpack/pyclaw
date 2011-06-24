"""
Module for PetClaw controller class.  The PetClaw controller is identical to the
PyClaw controller except for the default value of output_format.
"""

from pyclaw.controller import Controller as pyclawController

class Controller(pyclawController):
    def __init__(self):
        super(Controller,self).__init__()

        self.output_format = 'petsc'

    def sum_F_over_grid(self,state,i):
        return state.gFVec.strideNorm(i,0)

    def is_proc_0(self):
        from petsc4py import PETSc
        rank = PETSc.Comm.getRank(PETSc.COMM_WORLD)
        return rank == 0
