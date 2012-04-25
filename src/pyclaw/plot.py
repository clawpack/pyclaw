r"""Convenience routines for easily plotting with VisClaw."""

import os
import sys
import types

def plot(setplot_path=None,outdir="./_output",plotdir=None,htmlplot=False,iplot=False,
         file_format='ascii',**plot_kargs):
    r""""""
    
    # Construct a plot directory if not provided
    if plotdir is None:
        try: 
            plotdir = os.path.join(os.path.split(outdir)[:-2],"_plots")
        except AttributeError:
            plotdir = os.path.join(os.getcwd(),"_plots")
    
    if htmlplot or iplot:
        # Grab and import the setplot function
        # Use local setplot if available
        if setplot_path is None:
            local_setplot_path = os.path.join(os.getcwd(),'setplot.py')
            if os.path.exists(local_setplot_path):
                setplot_path = local_setplot_path

        # If setplot_path is still None then we import the default setplot
        if setplot_path is None:
            import visclaw.setplot_default as setplot_module
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
    plot(outdir=outdir,file_format=file_format,iplot=True)

def html_plot(outdir='./_output',file_format='ascii'):
    """
    Convenience function for creating html page with plots.
    """
    plot(outdir=outdir,file_format=file_format,htmlplot=True)
