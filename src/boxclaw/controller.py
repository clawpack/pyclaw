"""
Module for BoxClaw controller class.
"""

import fboxlib
from clawpack import pyclaw

class Controller(pyclaw.controller.Controller):
    """Parallel Controller Class"""

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.controller.Controller)

    def __init__(self):
        super(Controller,self).__init__()
        self.output_format = 'multifab'

    def is_proc_0(self):
        return fboxlib.mpi_rank() == 0

    def log_info(self, str):
        import logging
        if self.is_proc_0():
            logging.info(str)
        else:
            pass
