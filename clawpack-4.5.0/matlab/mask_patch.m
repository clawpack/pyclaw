function mask_patches(p,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh)

% Internal matlab routine for Clawpack graphics.

% Get stored Cartesian coordinates.  If we are not plotting a manifold or
% mapped grid, we probably don't need to store this, but we do it anyway
% so that we don't have always check for the manifold or mapped grid
% condition.

udata = get(p,'UserData');

% Use tolerance s so we don't mask cells that are just beyond masked patch.
% this value depends on the precision of dx,dy,dz stored in the fort.qXXXX
% files. If they are stored in single precision, than setting the
% following to 1e-8 will leave gaps in the masking.  So assuming single
% precision is probably a good idea here.

s = 1e-5;
if (s > min([udata.dx, udata.dy, udata.dz]))
  error('mask_patch : s is too big for dx, dy, dz\n');
end;

vc = udata.cartCoords;
[xdata_like,ydata_like, zdata_like]  = get_xyzlike(vc(:,1),vc(:,2), vc(:,3),sdir);
[xlow_like, ylow_like, zlow_like]    = get_xyzlike(xlow,ylow,zlow,sdir);
[xhigh_like, yhigh_like, zhigh_like] = get_xyzlike(xhigh,yhigh,zhigh,sdir);

nan_mask = ydata_like > ylow_like+s & ydata_like < yhigh_like-s & ...
    zdata_like > zlow_like+s & zdata_like < zhigh_like-s;

% Mask out physical vertices that should be hidden.
vp = get(p,'Vertices');
vp(nan_mask,1) = nan;
vp(nan_mask,2) = nan;
vp(nan_mask,3) = nan;

set(p,'Vertices',vp);
