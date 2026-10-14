from six.moves import range
def interactive_plot(outdir='./_output',file_format='petsc',setplot=None):
    """
    Convenience function for launching an interactive plotting session.
    """
    from clawpack.pyclaw.plot import plot
    plot(setplot,outdir=outdir,file_format=file_format,iplot=True,htmlplot=False)

def html_plot(outdir='./_output',file_format='petsc',setplot=None):
    """
    Convenience function for creating html page with plots.
    """
    from clawpack.pyclaw.plot import plot
    plot(setplot,outdir=outdir,file_format=file_format,htmlplot=True,iplot=False)

def plotPetsc(clawobj,delay=1):
    """
    Takes either a controller or solution object and prints each frame
    using PETSc.Viewer.
    """
    from petsc4py import PETSc
    from clawpack.pyclaw import controller, solution

    if isinstance(clawobj,controller.Controller):
        for n in range(0,clawobj.num_output_times):
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


