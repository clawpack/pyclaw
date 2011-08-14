#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing SharpClaw solvers for PetClaw

#  File:        sharpclaw.py
#  Created:     2010-03-20
#  Author:      David Ketcheson
"""
# ============================================================================
#      Copyright (C) 2010 David I. Ketcheson <david.ketcheson@kaust.edu.sa>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

from clawpack import PetSolver
from pyclaw.sharpclaw import SharpClawSolver1D, SharpClawSolver2D


class SharpClawSolver1D(PetSolver,SharpClawSolver1D):
    """
    
    1D parallel SharpClaw solver.

    Note that there are routines here for interfacing with the fortran time
    stepping routines and the python time stepping routines.  The ones used are
    dependent on the kernel_language argument given to the initialization of
    the solver (defaults to fortran).
    
    """
    def setup(self,solution):
        """
        Allocate RK stage arrays and fortran routine work arrays.
        """
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        state = solution.state
        state.set_stencil_width(self.mbc)
        # End hack

        self.allocate_rk_stages(solution)
        self.set_mthlim()
 
        if self.kernel_language=='Fortran':
            from sharpclaw1 import clawparams, workspace, reconstruct
            import sharpclaw1
            state = solution.states[0]
            state.set_cparam(sharpclaw1)
            self.set_fortran_parameters(state,clawparams,workspace,reconstruct)


class SharpClawSolver2D(PetSolver,SharpClawSolver2D):
    """
    
    2D parallel SharpClaw solver.  
    
    """

    def setup(self,solution):
        """
        Allocate RK stage arrays and fortran routine work arrays.
        """
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        state = solution.state
        state.set_stencil_width(self.mbc)
        # End hack

        self.allocate_rk_stages(solution)
        self.set_mthlim()
 
        if self.kernel_language=='Fortran':
            from sharpclaw2 import clawparams, workspace, reconstruct
            import sharpclaw2
            state.set_cparam(sharpclaw2)
            self.set_fortran_parameters(state,clawparams,workspace,reconstruct)


