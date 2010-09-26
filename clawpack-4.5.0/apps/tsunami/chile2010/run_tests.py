#! /usr/bin/python

"""
Run several tests with different parameters, in this case different numbers
refinement levels.
"""

from setrun import setrun
from setplot import setplot
from pyclaw.runclaw import runclaw
from pyclaw.plotters.plotclaw import plotclaw

# initialize rundata using setrun but then change some things for each run:
rundata = setrun()

#--------------
# Run 1:
#--------------
rundata.clawdata.mxnest = 1
rundata.write()

runclaw(xclawcmd = "xgeoclaw", outdir="_output_1level")
plotclaw(outdir="_output_1level", plotdir="_plots_1level")

#--------------
# Run 2:
#--------------
rundata.clawdata.mxnest = 2
rundata.geodata.wavetolerance = 0.01
rundata.write()

runclaw(xclawcmd = "xgeoclaw", outdir="_output_2level")
plotclaw(outdir="_output_2level", plotdir="_plots_2level")

#--------------
# Run 3:
#--------------
rundata.clawdata.mxnest = 3
rundata.geodata.wavetolerance = 0.1
rundata.write()

runclaw(xclawcmd = "xgeoclaw", outdir="_output_3level")
plotclaw(outdir="_output_3level", plotdir="_plots_3level")

print "Plots should be in _plots_1level, _plots_2level, _plots_3level"

