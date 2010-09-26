function [isect,qcm2] = interp_data_3d(xc,yc,zc,xe,ye,ze,qcm3,...
    sdir,sval,interpmethod)

% Internal matlab routine for Clawpack graphics.

% --------------------------------------
% First Check to see that 3D patch intersects slice n.
% --------------------------------------

vlow  = [xe(1), ye(1), ze(1)];
vhigh = [xe(end),ye(end),ze(end)];

s = 1e-6;
idir = findstr(lower(sdir),'xyz');
isect = (vlow(idir)-s <= sval) & (sval <= vhigh(idir)+s);
if (isect == 0)
  % This slice doesn't intersect data.
  qcm2 = [];
  return;
end;

% --------------------------------------
% Now setup data for interpolation
% --------------------------------------

[xc_like,yc_like, zc_like] = get_xyzlike(xc,yc,zc,sdir);
[xe_like,ye_like, ze_like] = get_xyzlike(xe,ye,ze,sdir);

% Make sure slice constant isn't outside of (x,y,z) data range, or we'll get
% NANs when we really want number.
sval_interp = min([max([sval,xc_like(1)]),xc_like(end)]);

% Set up 2d data for interpolation.
[ycm2_like,zcm2_like] = meshgrid(yc_like,zc_like);
xcm2_like = 0*ycm2_like + sval_interp;

% Recover phsical (x,y,z) coordinates for intepolation.
[xcm2, ycm2, zcm2] = get_xyz(xcm2_like, ycm2_like, zcm2_like,sdir);

% Get 3d cubes for centers of each cell.
[xcm3,ycm3,zcm3] = meshgrid(xc,yc,zc);

qcm2 = interp3(xcm3,ycm3,zcm3,qcm3,xcm2,ycm2,zcm2,interpmethod);
