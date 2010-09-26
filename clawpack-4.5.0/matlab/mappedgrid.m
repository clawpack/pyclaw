% MappedGrid parameter for plotting using a non-uniform grid.
%
%     If MappedGrid = 1 in SETPLOT1, SETPLOT2 or SETPLOT3 matlab scripts,
%     the Clawpack graphics routines will call a user-defined script
%     'mapc2p.m', which specifies mapping from Cartesian coordinates to
%     physical coordinates.  The script has the following input and output
%     arguments :
%
%     In 1d :
%           xp = mapc2p(xc);
%
%     In 2d :
%           [xp,yp] = mapc2p(xc,yc);
%
%     In 3d :
%           [xp,yp,zp] = mapc2p(xc,yc,zc);
%
%     Example :
%
%           function [xp,yp] = mapc2p(xc,yc);
%           % For polar plots in 2d.
%           xp = yc*cos(2*pi*xc);
%           yp = yc*sin(2*pi*xc);
%
%     This script should be in the current directory, or in the Matlab path.
%
%     When plotting results of a multiblock computation, the user may want
%     to use the GETBLOCKNUMBER function to specify how the mapping should
%     behave for each block.
%
%     For Scatter or Line plots in 2d or 3d (PlotType = 4), the MappedGrid
%     option will be used in the following manner :  If the user has set
%     (x0,y0,z0), then the distance to these points will be computed using
%     physical coordinates, e.g. in 3d
%
%          [xp,yp,zp] = mapc2p(xc,yc,zc);
%          r = sqrt((xp-x0).^2 + (yp-y0).^2 +  (zp-z0).^2);
%
%     Input arguments to a user-defined 'map1d' function are in Cartesian
%     coordinates. The user should call 'mapc2p' from within their 'map1d'
%     routine.
%
%     See  also SETPLOT, MANIFOLD, GETBLOCKNUMBER.

error(['mappedgrid : This is not a stand-alone routine.  Set MappedGrid in ',...
  'setplot1, 2, or 3 to use mapped grids\n']);
