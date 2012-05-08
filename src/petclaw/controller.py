"""
Module for PetClaw controller class.  The PetClaw controller is identical to the
PyClaw controller except for the default value of output_format.
"""

from clawpack.pyclaw.controller import Controller as pyclawController

class Controller(pyclawController):
    def __init__(self):
        super(Controller,self).__init__()

        self.output_format = 'petsc'

    def is_proc_0(self):
        from petsc4py import PETSc
        rank = PETSc.Comm.getRank(PETSc.COMM_WORLD)
        return rank == 0

    def log_info(self, str):
        import logging
        if self.is_proc_0():
            logging.info(str)
        else:
            pass
