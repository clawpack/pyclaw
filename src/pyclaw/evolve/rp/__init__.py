#!/usr/bin/env python
# encoding: utf-8
"""
__init__.py

Created by Kyle Mandli on 2008-09-15.
Copyright (c) 2008 University of Washington. All rights reserved.
"""

rp_solver_list_1d = ['advection','vc_advection','acoustics','burgers','shallow_roe',
                  'shallow_hll','shallow_exact','euler_roe','nel']
rp_solver_list_2d = []
rp_solver_list_3d = []

# Import 1d Riemann solvers
from rp_advection import rp_advection_1d
from rp_vc_advection import rp_vc_advection_1d
from rp_acoustics import rp_acoustics_1d
from rp_burgers import rp_burgers_1d
from rp_shallow import rp_shallow_roe_1d, rp_shallow_hll_1d, rp_shallow_exact_1d
from rp_euler import rp_euler_roe_1d
from rp_nel import rp_nel_1d
