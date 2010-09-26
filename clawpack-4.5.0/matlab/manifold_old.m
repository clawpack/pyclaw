% Manifold parameter for plotting 2d data on a manifold.
%
%    If Manifold = 1 in the SETPLOT2 matlab scripts, the Clawpack graphics
%    routines will call a user-defined script 'mapc2m.m', which specifies
%    mapping from 2d Cartesian coordinates to 3d physical coordinates.
%    The script has the following input and output arguments :
%
%           [xp,yp,zp] = mapc2m(xc,yc);
%
%     Example :
%
%           function [xp,yp,zp] = mapc2m(xc,yc);
%           % For plotting data on the surface of a unit sphere.
%           r = 1;
%           xp = r*cos(yc).*cos(xc);
%           yp = r*cos(yc).*sin(xc);
%           zp = r*sin(yc);
%
%     This script should be in the current directory, or in the Matlab path.
%
%     See  also SETPLOT, MAPPEDGRID.
%

error(['manifold : This routine cannot be called.  Set Manifold in ',...
  'setplot2 to plot data on a manifold.\n']);
