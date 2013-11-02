r"""Convenience routines for easily plotting with VisClaw."""

import os
import sys
import types

def plot(setplot=None,outdir="./_output",plotdir=None,htmlplot=False,iplot=True,
         file_format='ascii',**plot_kargs):
    r"""
    setplot can be a function or a path to a file.
    """
    
    # Construct a plot directory if not provided
    if plotdir is None:
        try: 
            plotdir = os.path.join(os.path.split(outdir)[:-2],"_plots")
        except AttributeError:
            plotdir = os.path.join(os.getcwd(),"_plots")
    
    if htmlplot or iplot:
        # Grab and import the setplot function
        if not isinstance(setplot,types.FunctionType):
            # Use local setplot if available
            if setplot is None:
                local_setplot_path = os.path.join(os.getcwd(),'setplot.py')
                if os.path.exists(local_setplot_path):
                    setplot = local_setplot_path

            # If setplot is still None then we import the default setplot
            if setplot is None:
                import clawpack.visclaw.setplot_default as setplot_module
            else:
                path = os.path.abspath(os.path.expandvars(os.path.expanduser(setplot)))
                setplot_module_dir = os.path.dirname(path)
                setplot_module_name = os.path.splitext(os.path.basename(setplot))[0]
                sys.path.insert(0,setplot_module_dir)
                setplot_module = __import__(setplot_module_name)
            reload(setplot_module)
            setplot = lambda plotdata:setplot_module.setplot(plotdata,**plot_kargs)
        
        if not isinstance(setplot,types.FunctionType):
            raise ImportError("Failed to import %s.setplot" % setplot_module_name)
        
        if iplot:
            from clawpack.visclaw import Iplotclaw
        
            ip = Iplotclaw.Iplotclaw(setplot=setplot,outdir=outdir)
            ip.plotdata.format = file_format
        
            ip.plotloop()
            
        if htmlplot:
            from clawpack.visclaw import plotclaw
            plotclaw.plotclaw(outdir,plotdir,format=file_format,setplot=setplot)
        

# These now just point to the above more generic function
def interactive_plot(outdir='./_output',file_format='ascii',setplot=None):
    """
    Convenience function for launching an interactive plotting session.
    """
    plot(setplot,outdir=outdir,file_format=file_format,iplot=True,htmlplot=False)

def html_plot(outdir='./_output',file_format='ascii',setplot=None):
    """
    Convenience function for creating html page with plots.
    """
    plot(setplot,outdir=outdir,file_format=file_format,htmlplot=True,iplot=False)
