r"""Convenience routines for easily plotting with VisClaw."""

def interactive_plot(outdir='./_output',file_format='ascii'):
    """
    Convenience function for launching an interactive plotting session.
    """
    from visclaw import Iplotclaw
    ip=Iplotclaw.Iplotclaw()
    ip.plotdata.outdir=outdir
    ip.plotdata.format=file_format
    ip.plotloop()

def html_plot(outdir='./_output',file_format='ascii'):
    """
    Convenience function for creating html page with plots.
    """
    from visclaw import plotclaw
    plotclaw.plotclaw(outdir,format=file_format)
