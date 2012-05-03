def interactive_plot(outdir='./_output',format='petsc'):
    """
    Convenience function for launching an interactive plotting session.
    """
    import clawpack.visclaw.Iplotclaw as Iplotclaw
    ip=Iplotclaw.Iplotclaw()
    ip.plotdata.outdir=outdir
    ip.plotdata.format=format
    ip.plotloop()

def html_plot(outdir='./_output',format='petsc'):
    """
    Convenience function for creating html page with plots.
    """
    import clawpack.visclaw.plotclaw as plotclaw
    plotclaw.plotclaw(outdir,format=format)

def plotPetsc(clawobj,delay=1):
    """
    Takes either a controller or solution object and prints each frame
    using PETSc.Viewer.
    """
    from petsc4py import PETSc
    from clawpack.pyclaw import controller, solution

    if isinstance(clawobj,controller.Controller):
        for n in xrange(0,clawobj.num_output_times):
            sol = clawobj.frames[n]
            viewer = PETSc.Viewer.DRAW(sol.patch.gqVec.comm)
            OptDB = PETSc.Options()
            OptDB['draw_pause'] = delay
            viewer(sol.patch.gqVec)

    elif isinstance(clawobj,solution.Solution):
        viewer = PETSc.Viewer.DRAW(clawobj.patch.gqVec.comm)
        OptDB = PETSc.Options()
        OptDB['draw_pause'] = -1
        viewer(clawobj.patch.gqVec)


