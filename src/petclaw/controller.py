"""
Module for PetClaw controller class.
"""
from clawpack import pyclaw

class Controller(pyclaw.controller.Controller):
    """ Parallel Controller Class

    Defaults to petsc output_format, logs only from process 0.
    """
    
    __doc__ += pyclaw.util.add_parent_doc(pyclaw.controller.Controller)

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
