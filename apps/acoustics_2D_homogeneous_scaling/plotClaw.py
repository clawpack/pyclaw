


from pyclaw.plotters import Iplotclaw
ip=Iplotclaw.Iplotclaw()
ip.plotdata.outdir='./_output/'
ip.plotdata.format='petsc'
ip.plotloop()
