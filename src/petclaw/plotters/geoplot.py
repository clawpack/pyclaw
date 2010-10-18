"""
Useful things for plotting GeoClaw results.
"""

from petclaw.plotters import colormaps

# Colormaps from geoclaw
# Color attributes, single instance per run
# Colors
black = [0.0,0.0,0.0]
white = [1.0,1.0,1.0]
red = [1.0,0.0,0.0]
green = [0.0,1.0,0.0];
dark_green = [0.1,0.4,0.0];
light_green = [0.8,1.0,0.5];
blue = [0.0,0.0,1.0];
dark_blue = [0.2,0.2,0.7];
light_blue = [0.5,0.5,1.0];
blue_green = [0.0,1.0,1.0];
tan = [0.9,0.8,0.2];
tan = [0.8,0.5,0.2];
brown = [0.9,0.8,0.2];
gray8 = [0.8,0.8,0.8];
purple = [0.8,0.3,0.8];

# Colormaps
TSUNAMI_MAX_AMPLITUDE = 0.6
tsunami_colormap = colormaps.make_colormap({-TSUNAMI_MAX_AMPLITUDE:blue,
                                            0.0:blue_green,
                                            TSUNAMI_MAX_AMPLITUDE:red})
                                            
land1_colormap = colormaps.make_colormap({0.0:dark_green,
                                          1000.0:green,
                                          2000.0:light_green,
                                          4000.0:tan})
                                         
land2_colormap = colormaps.make_colormap({0:dark_green,
                                          50:green,
                                          100:light_green,
                                          200:tan})
                                          
water_land_colormap = colormaps.make_colormap({-1000:dark_blue,
                                               -500:blue,
                                               0:light_blue,
                                               .1:tan,
                                               5:tan,
                                               6:dark_green,
                                               1000:green,
                                               2000:light_green,
                                               4000:tan})
                                               
bathy1_colormap = colormaps.make_colormap({-1000:brown,
                                           0:tan,
                                           .1:dark_green,
                                           1000:green,
                                           2000:light_green})
       
bathy2_colormap = colormaps.make_colormap({-1000:brown,
                                           -100:tan,
                                           0:dark_green,
                                           .1:dark_green,
                                           1000:green,
                                           2000:light_green})

bathy3_colormap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
                                           -0.01:[0.95,0.9,0.7],
                                           .01:[.5,.7,0],
                                           1:[.2,.5,.2]})

                                           
colormaps_list = {"tsunami":tsunami_colormap,"land1":land1_colormap,
             "land2":land2_colormap,"water_land":water_land_colormap,
             "bathy1":bathy1_colormap,"bathy2":bathy2_colormap}
        
def plot_colormaps():
    r"""Plots all colormaps avaiable or the ones specified"""
        
    import numpy as np
    import matplotlib.pyplot as plt
    
    a = np.linspace(0, 1, 256).reshape(1,-1)
    a = np.vstack((a,a))
    
    nmaps = len(colormaps_list) + 1

    fig = plt.figure(figsize=(5,10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    
    for i,name in enumerate(colormaps_list):
        ax = plt.subplot(nmaps,1,i+1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=colormaps_list[name], origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10, horizontalalignment='right')

    #plt.show()

land_colors = colormaps.make_colormap({0:[.5,.7,0], 1:[.2,.5,.2]})
# water_colors = colormaps.make_colormap({-1.:'r', 0.:[0, .8, .8], 1.: 'b'})
# land_colors = land2_colormap
water_colors = tsunami_colormap

# Plotting functions

# The drytol parameter is used in masking land and water and
# affects what color map is used for cells with small water depth h.
# The best value to use often depends on the application and can
# be set for an application by setting current_data.user.drytol in
# a beforeframe function, for example.  If it's not set by the user,
# the following default value is used (in meters):

drytol_default = 1.e-3


def topo(current_data):
   """ 
   Return topography = eta - h. 
   Surface eta is assumed to be output as 4th column of fort.q files.
   """
   q = current_data.q
   h = q[:,:,0]
   eta = q[:,:,3]
   topo = eta - h
   return topo

def land(current_data):
   """
   Return a masked array containing the surface elevation only in dry cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[:,:,0]
   eta = q[:,:,3]
   land = ma.masked_where(h>drytol, eta)
   return land

def water(current_data):
   """Deprecated: use surface instead."""
   from numpy import ma
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[:,:,0]
   eta = q[:,:,3]
   water = ma.masked_where(h<=drytol, eta)
   return water

def depth(current_data):
   """
   Return a masked array containing the depth of fluid only in wet cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[:,:,0]
   depth = ma.masked_where(h<=drytol, h)
   return depth

def surface(current_data):
   """
   Return a masked array containing the surface elevation only in wet cells.
   Surface is eta = h+topo, assumed to be output as 4th column of fort.q
   files.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[:,:,0]
   eta = q[:,:,3]
   water = ma.masked_where(h<=drytol, eta)
   return water

def surface_or_depth(current_data):
   """
   Return a masked array containing the surface elevation where the topo is 
   below sea level or the water depth where the topo is above sea level.
   Mask out dry cells.  Assumes sea level is at topo=0.
   Surface is eta = h+topo, assumed to be output as 4th column of fort.q
   files.
   """
   from numpy import ma, where
   drytol = getattr(current_data.user, 'drytol', drytol_default)
   q = current_data.q
   h = q[:,:,0]
   eta = q[:,:,3]
   topo = eta - h
   surface = ma.masked_where(h<=drytol, eta)
   depth = ma.masked_where(h<=drytol, h)
   surface_or_depth = where(topo<0, surface, depth)
   return surface_or_depth


def plot_topo_file(topoplotdata):
    """
    Read in a topo or bathy file and produce a pcolor map.
    """

    import os
    import pylab
    from petclaw.data import Data

    fname = topoplotdata.fname 
    topotype = getattr(topoplotdata, 'topotype', 3)
    cmap = getattr(topoplotdata, 'cmap', bathy3_colormap)
    climits = getattr(topoplotdata, 'climits', [-50,50])
    figno = getattr(topoplotdata, 'figno', 200)
    addcolorbar = getattr(topoplotdata, 'addcolorbar', True)
    addcontour = getattr(topoplotdata, 'addcontour', False)
    contour_levels = getattr(topoplotdata, 'contour_levels', [0, 0])
    xlimits = getattr(topoplotdata, 'xlimits', None)
    ylimits = getattr(topoplotdata, 'ylimits', None)
    coarsen = getattr(topoplotdata, 'coarsen', 1)

    if abs(topotype) != 3:
        print "*** Only topotype = 3, -3 implemented so far."
        return

    file = open(fname, 'r')
    lines = file.readlines()
    ncols = int(lines[0].split()[0])
    nrows = int(lines[1].split()[0])
    xllcorner = float(lines[2].split()[0])
    yllcorner = float(lines[3].split()[0])
    cellsize = float(lines[4].split()[0])
    NODATA_value = int(lines[5].split()[0])

    print "Loading file ",fname
    print "   nrows = %i, ncols = %i" % (nrows,ncols)
    topo = pylab.loadtxt(fname,skiprows=6,dtype=float)
    print "   Done loading"

    if 0:
        topo = []
        for i in range(nrows):
            topo.append(pylab.array(lines[6+i],))
        print '+++ topo = ',topo
        topo = pylab.array(topo)

    topo = pylab.flipud(topo)
    if topotype < 0:
        topo = -topo

    x = pylab.linspace(xllcorner, xllcorner+ncols*cellsize, ncols)
    y = pylab.linspace(yllcorner, yllcorner+nrows*cellsize, nrows)
    print "Shape of x, y, topo: ", x.shape, y.shape, topo.shape

    if coarsen > 1:
        topo = topo[slice(0,nrows,coarsen), slice(0,ncols,coarsen)]
        x = x[slice(0,ncols,coarsen)]
        y = y[slice(0,nrows,coarsen)]
        print "Shapes after coarsening: ", x.shape, y.shape, topo.shape


    #import pdb; pdb.set_trace()

    if figno:
        pylab.figure(figno)

    pylab.pcolor(x,y,topo,cmap=cmap)
    pylab.clim(climits)
    pylab.axis('scaled')

    if addcolorbar:
        pylab.colorbar()

    if addcontour:
        pylab.contour(x,y,topo,levels=contour_levels,colors='k')

    gridedges_show = True
    if gridedges_show:
        pylab.plot([x[0],x[-1]],[y[0],y[0]],'k')
        pylab.plot([x[0],x[-1]],[y[-1],y[-1]],'k')
        pylab.plot([x[0],x[0]],[y[0],y[-1]],'k')
        pylab.plot([x[-1],x[-1]],[y[0],y[-1]],'k')

    fname2 = os.path.splitext(fname)[0]
    pylab.text(xllcorner+cellsize, yllcorner+cellsize, fname2, color='m')

    topodata = Data()
    topodata.x = x
    topodata.y = y
    topodata.topo = topo

    return topodata
