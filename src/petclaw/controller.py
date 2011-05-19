"""
Module for PetClaw controller class.  The PetClaw controller is identical to the
PyClaw controller except for the default value of output_format.
"""

from pyclaw.controller import Controller as pyclawController

class Controller(pyclawController):
    def __init__(self):
        super(Controller,self).__init__()

        print self.xdir
        self.output_format = 'petsc'
