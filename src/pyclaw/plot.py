r"""Convenience routines for easily plotting with VisClaw."""

def plot(setplot_path=None,outdir="./_output",plotdir=None,htmlplot=False,iplot=False,
         file_format='ascii',**plot_kargs):
    r""""""
    
    # Construct a plot directory if not provided
    if plotdir is None:
        plotdir = os.path.join(os.path.split(outdir)[:-2],"_plots")
    
    if htmlplot or iplot:
        # Grab and import the setplot function
        if setplot_path is None:
            setplot_module = __import__("visclaw.setplot_default")
        else:
            path = os.path.abspath(os.path.expandvars(os.path.expanduser(setplot_path)))
            setplot_module_dir = os.path.dirname(path)
            setplot_module_name = os.path.splitext(os.path.basename(setplot_path))[0]
            sys.path.insert(0,setplot_module_dir)
        setplot_module = __import__(setplot_module_name)
        reload(setplot_module)
        setplot = lambda plotdata:setplot_module.setplot(plotdata,**plot_kargs)
        
        if not isinstance(setplot,types.FunctionType):
            raise ImportError("Failed to import %s.setplot" % setplot_module_name)
        
        if iplot:
            from visclaw import Iplotclaw
        
            ip = Iplotclaw.Iplotclaw(setplot=setplot)
            ip.plotdata.outdir = outdir
            ip.plotdata.format = file_format
        
            ip.plotloop()
            
        if htmlplot:
            from visclaw import plotclaw            
            plotclaw.plotclaw(outdir,plotdir,format=file_format,setplot=setplot)
        

# These now just point to the above more generic function
def interactive_plot(outdir='./_output',file_format='ascii'):
    """
    Convenience function for launching an interactive plotting session.
    """
    plot(outdir=outdir,file_format=file_format)

def html_plot(outdir='./_output',file_format='ascii'):
    """
    Convenience function for creating html page with plots.
    """
    plot(outdir=outdir,file_format=file_format)
