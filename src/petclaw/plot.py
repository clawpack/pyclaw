def plotPetsc(clawobj,delay=1):
    """
    Takes either a controller or solution object and prints each frame
    using PETSc.Viewer.
    """
    from petsc4py import PETSc
    import pyclaw.controller, pyclaw.solution

    if isinstance(clawobj,pyclaw.controller.Controller):
        for n in xrange(0,clawobj.nout):
            sol = clawobj.frames[n]
            plotTitle="time: {0}".format(sol.t)
            viewer = PETSc.Viewer.DRAW(sol.grid.gqVec.comm)
            OptDB = PETSc.Options()
            OptDB['draw_pause'] = delay
            viewer(sol.grid.gqVec)

    elif isinstance(clawobj,pyclaw.solution.Solutuion):
        viewer = PETSc.Viewer.DRAW(clawobj.grid.gqVec.comm)
        OptDB = PETSc.Options()
        OptDB['draw_pause'] = -1
        viewer(clawobj.grid.gqVec)

def plotInteractive(outdir='./_output'):
    """
    Convenience function for launching an interactive plotting session.
    """
    from pyclaw.plotters import Iplotclaw
    ip=Iplotclaw.Iplotclaw()
    ip.plotdata.outdir=outdir
    ip.plotdata.format='petsc'
    ip.plotloop()

def plotHTML(outdir='./_output'):
    """
    Convenience function for creating html page with plots.
    """
    from pyclaw.plotters import plotclaw
    plotclaw.plotclaw('./_output',format='petsc')
