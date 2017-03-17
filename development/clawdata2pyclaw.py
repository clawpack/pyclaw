"""
Script to convert Clawpack 4.6.x problem setup to PyClaw setup.
Automatically writes a basic PyClaw script with all the options
from the setrun.py in the current directory.  Also generates a 
Makefile and a wrapper for qinit.f.  Additional wrappers will be
needed for any other custom Fortan code, such as setaux, b4step, etc.

If you try this out, please raise issues in the PyClaw tracker for 
anything that doesn't work.
"""
from __future__ import absolute_import
from __future__ import print_function
import setrun
import six

rundata=setrun.setrun()
clawdata=rundata.clawdata
probdata=rundata.probdata
ndim = clawdata.ndim

print("Writing run.py...")
outfile = 'pyclaw/run.py'
output=open(outfile,'w')

output.write("#!/usr/bin/env python\n")
output.write("# encoding: utf-8\n\n")


# write wrappers for Fortran functions
output.write("def fortran_qinit_wrapper(solver,state):\n")
output.write('    """\nWraps Fortran routine qinit.f"""\n')
output.write("    grid = state.grid\n")
output.write("    meqn = state.num_eqn\n")
output.write("    maux = state.num_aux\n")
output.write("    mbc = solver.num_ghost\n")
output.write("    q = state.q\n")
output.write("    aux = state.aux\n")
output.write("    t = state.t\n")
output.write("    mx = grid.num_cells[0]\n")
output.write("    dx = grid.delta[0]\n")
output.write("    xlower = grid.lower[0]\n")
if ndim>1:
    output.write("    my = grid.num_cells[1]\n")
    output.write("    dy = grid.delta[1]\n")
    output.write("    ylower = grid.lower[1]\n")
if ndim>2:
    output.write("    mz = grid.num_cells[2]\n")
    output.write("    dz = grid.delta[2]\n")
    output.write("    zlower = grid.lower[2]\n")

output.write("\n")
output.write("    import problem\n")
if ndim==1:
    output.write("    state.q = problem.qinit(mx,meqn,mbc,mx,xlower,dx,maux,aux)")
elif ndim==2:
    output.write("    state.q = problem.qinit(mx,meqn,mbc,mx,my,xlower,ylower,dx,dy,maux,aux)")
elif ndim==3:
    output.write("    state.q = problem.qinit(mx,meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,maux,aux)")
output.write("\n\n\n")


output.write("# Start main script")
output.write("import numpy as np\n")
output.write("import pyclaw\n\n")

output.write("solver = pyclaw.ClawSolver%sD()\n" % ndim)
output.write("import riemann\n")
output.write("solver.rp = riemann.<RIEMANN SOLVER NAME HERE>\n\n")

import pyclaw
exec("solver = pyclaw.ClawSolver%sD()" % ndim)
output.write("#Set all solver attributes\n")
solver_attrs = clawdata.__dict__
for key,value in six.iteritems(solver_attrs):
    if hasattr(solver,key):
        output.write("solver.%s = %s\n" % (key,value))
output.write("\n")

output.write("""\n\
# Choice of BCs at lower and upper:\n\
#   0 => user specified (must modify bcN.f to use this option)\n\
#   1 => extrapolation (non-reflecting outflow)\n\
#   2 => periodic (must specify this at both boundaries)\n\
#   3 => solid wall for systems where q(2) is normal velocity\n""")

output.write("solver.bc_lower[0] = %s\n" % clawdata.bc_xlower)
output.write("solver.bc_upper[0] = %s\n" % clawdata.bc_xupper)
if ndim>1:
    output.write("solver.bc_lower[1] = %s\n" % clawdata.bc_ylower)
    output.write("solver.bc_upper[1] = %s\n" % clawdata.bc_yupper)
if ndim>2:
    output.write("solver.bc_lower[2] = %s\n" % clawdata.bc_zlower)
    output.write("solver.bc_upper[2] = %s\n" % clawdata.bc_zupper)
output.write("\n")

output.write("solver.max_steps = %s\n" % clawdata.steps_max)
output.write("solver.num_waves = %s\n" % clawdata.mwaves)
output.write("solver.limiters = %s\n" % clawdata.limiter)
output.write("\n")

output.write("""
# Source terms splitting:\n\
#   src_split == 0  => no source term (src routine never called)\n\
#   src_split == 1  => Godunov (1st order) splitting used, \n\
#   src_split == 2  => Strang (2nd order) splitting used,  not recommended.\n""")
output.write("solver.source_split = %s\n" % clawdata.src_split)
output.write("\n")

output.write("# Initialize domain\n")
output.write("x = pyclaw.Dimension(%s,%s,%s)\n" % (clawdata.xlower,clawdata.xupper,clawdata.mx))
if ndim>1:
    output.write("y = pyclaw.Dimension(%s,%s,%s)\n" % (clawdata.ylower,clawdata.yupper,clawdata.my))
if ndim>2:
    output.write("z = pyclaw.Dimension(%s,%s,%s)\n" % (clawdata.zlower,clawdata.zupper,clawdata.mz))

output.write("domain = pyclaw.Domain")
if ndim == 1:
    output.write("(x)")
elif ndim == 2:
    output.write("(x,y)")
elif ndim == 3:
    output.write("(x,y,z)")
output.write("\n\n")

output.write("# Initialize state\n")
output.write("num_eqn = %s\n" % clawdata.meqn)
output.write("num_aux = %s\n" % clawdata.maux)
output.write("state = pyclaw.State(domain,num_eqn,num_aux)\n")
output.write("state.capa_index = %s" % clawdata.mcapa)

output.write("# Set problem data\n")
for param,value in six.iteritems(probdata):
    output.write("state.problem_data['%s']=%s\n" % (param,value))

output.write("\n")

output.write("# Initialize controller and solution\n")
output.write("claw = pyclaw.Controller()\n")
output.write("claw.solution = pyclaw.Solution(state,domain)\n")
output.write("claw.solution.t = %s\n" % clawdata.t0)
output.write("claw.solver = solver\n")
output.write("claw.tfinal = %s\n" % clawdata.output_tfinal)
output.write("claw.output_style = %s\n" % clawdata.output_style)
output.write("claw.num_output_times = %s\n" % clawdata.output_ntimes)
output.write("claw.verbosity = %s\n" % clawdata.verbosity)

output.write("\nstatus = claw.run()\n\n\n")

output.close()


# Write Makefile
print("Writing Makefile...")
makefile = 'pyclaw/Makefile' # To avoid stomping existing Makefile
outmake=open(makefile,'w')

outmake.write("all:\n\tmake classic1.so\n\tmake problem.so\n\n")

outmake.write("# Put all problem-specific Fortran files here:\n")
outmake.write("problem.so: qinit.f setaux.f src.f mapc2p.f\n\t")
outmake.write(r"""$(F2PY) -m problem -c $^""")
outmake.write("\n\n")
outmake.write(r"""include $(PYCLAW)/Makefile.common""")
outmake.write("\n\n")
outmake.close()

print("""Don't forget to do the following manually:\n\
         1. interleave your Fortran code \n\
         2. Fill in the Riemann solver in run.py\n\
         3. Add 'cf2py intent(out) q' to your qinit.f\n\
         3. write any code needed to replicate setprob.f functionality.""")
