def interactive_plot(outdir='./_output',file_format='ascii',setplot=None):
    """
    Convenience function for launching an interactive plotting session.
    """
    from clawpack.pyclaw.plot import plot
    plot(setplot,outdir=outdir,file_format=file_format,iplot=True,htmlplot=False)

def html_plot(outdir='./_output',file_format='ascii',setplot=None):
    """
    Convenience function for creating html page with plots.
    """
    from clawpack.pyclaw.plot import plot
    plot(setplot,outdir=outdir,file_format=file_format,htmlplot=True,iplot=False)
