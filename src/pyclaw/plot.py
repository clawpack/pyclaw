r"""Convenience routines for easily plotting with VisClaw."""

from __future__ import absolute_import
import os
import sys
import types
import six

def plot(setplot=None, outdir="./_output", plotdir=None, htmlplot=False, 
         iplot=True, file_format='ascii', **plot_kargs):
    r"""setplot can be a function or a path to a file."""
    
    # Construct a plot directory if not provided
    if plotdir is None:
        try: 
            plotdir = os.path.join(outdir,"../_plots")
        except AttributeError:
            plotdir = os.path.join(os.getcwd(),"_plots")
    
    if htmlplot or iplot:
        setplot_func = None
        # No setplot specified, try to use a local file
        if setplot is None:
            # Grab and import the setplot function
            local_setplot_path = os.path.join(os.getcwd(),'setplot.py')
            if os.path.exists(local_setplot_path):
                setplot = local_setplot_path

        # Fetch setplot function depending on type of setplot
        if isinstance(setplot, types.FunctionType):
            # setplot points to a function
            setplot_func = lambda plotdata:setplot(plotdata, **plot_kargs)

        elif isinstance(setplot, types.ModuleType):
            # setplot points to a module
            setplot_func = lambda plotdata:setplot.setplot(plotdata, **plot_kargs)
            
        elif isinstance(setplot, six.string_types):
            # setplot contains a path to a module
            path = os.path.abspath(os.path.expandvars(os.path.expanduser(setplot)))
            setplot_module_dir = os.path.dirname(path)
            setplot_module_name = os.path.splitext(os.path.basename(setplot))[0]
            sys.path.insert(0,setplot_module_dir)
            setplot_module = __import__(setplot_module_name)
            setplot_func = lambda plotdata:setplot_module.setplot(plotdata, **plot_kargs)
        
        if not isinstance(setplot_func, types.FunctionType):
            # Everything else has failed, use default setplot
            import clawpack.visclaw.setplot_default as setplot_module
            setplot_func = setplot_module.setplot

        # Interactive plotting
        if iplot:
            from clawpack.visclaw import Iplotclaw
        
            ip = Iplotclaw.Iplotclaw(setplot=setplot_func, outdir=outdir)
            ip.plotdata.format = file_format
        
            ip.plotloop()
            
        # Static HTML plotting
        if htmlplot:
            from clawpack.visclaw import plotclaw
            plotclaw.plotclaw(outdir, plotdir, format=file_format,
                                               setplot=setplot_func)
        

# These now just point to the above more generic function
def interactive_plot(outdir='./_output', file_format='ascii', setplot=None):
    """Convenience function for launching an interactive plotting session."""
    plot(setplot, outdir=outdir, file_format=file_format, iplot=True, 
                  htmlplot=False)


def html_plot(outdir='./_output', file_format='ascii', setplot=None):
    """Convenience function for creating html page with plots."""
    plot(setplot, outdir=outdir, file_format=file_format, htmlplot=True, 
                  iplot=False)
