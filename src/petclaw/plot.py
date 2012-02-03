def interactive_plot(outdir='./_output',format='petsc'):
    """
    Convenience function for launching an interactive plotting session.
    """
    from visclaw import Iplotclaw
    ip=Iplotclaw.Iplotclaw()
    ip.plotdata.outdir=outdir
    ip.plotdata.format=format
    ip.plotloop()

def html_plot(outdir='./_output',format='petsc'):
    """
    Convenience function for creating html page with plots.
    """
    from visclaw import plotclaw
    plotclaw.plotclaw(outdir,format=format)

def plotPetsc(clawobj,delay=1):
    """
    Takes either a controller or solution object and prints each frame
    using PETSc.Viewer.
    """
    from petsc4py import PETSc
    import pyclaw.controller, pyclaw.solution

    if isinstance(clawobj,pyclaw.controller.Controller):
        for n in xrange(0,clawobj.num_output_times):
            sol = clawobj.frames[n]
            viewer = PETSc.Viewer.DRAW(sol.patch.gqVec.comm)
            OptDB = PETSc.Options()
            OptDB['draw_pause'] = delay
            viewer(sol.patch.gqVec)

    elif isinstance(clawobj,pyclaw.solution.Solution):
        viewer = PETSc.Viewer.DRAW(clawobj.patch.gqVec.comm)
        OptDB = PETSc.Options()
        OptDB['draw_pause'] = -1
        viewer(clawobj.patch.gqVec)


