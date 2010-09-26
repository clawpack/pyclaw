% From one-dimensional vectors xcenter and ycenter defining a rectangular grid
% in computational space, make 2d arrays xp, yp so that (xp(i,j), yp(i,j))
% is the corresponding point in physical space.  
%
% Assumes plotclaw2 has already been started so xcenter and ycenter are set.
%

[xc,yc] = meshgrid(xcenter,ycenter);

if MappedGrid
   [xp,yp] = mapc2p(xc,yc);
   end
