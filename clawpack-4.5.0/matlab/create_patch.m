function p = create_patch(xc,yc,zc,xe,ye,ze,qcm2,sdir,sval,...
    contourlevels,mappedgrid,manifold,grid_number)

% Internal matlab routine for Clawpack graphics.

[xe_like, ye_like, ze_like] = get_xyzlike(xe,ye,ze,sdir);
[xc_like, yc_like, zc_like] = get_xyzlike(xc,yc,zc,sdir);

% -----------x---------------------------------------
% Create patch with q data
% --------------------------------------------------
[yem2_like,zem2_like] = meshgrid(ye_like,ze_like);
xem2_like = 0*yem2_like + sval;
[xem2, yem2, zem2] = get_xyz(xem2_like,yem2_like, zem2_like,sdir);

% This command needs phsical x,y,z coordinates, not permuted coordinates.
p = patch(surf2patch(xem2,yem2,zem2));

% Now set up colors...
[ycm2_like,zcm2_like] = meshgrid(yc_like,zc_like);
xcm2_like = 0*ycm2_like + sval;
[xcm2, ycm2, zcm2] = get_xyz(xcm2_like, ycm2_like, zcm2_like, sdir);

% This needs physical coordinates as well.
setcolors(p,xcm2, ycm2, zcm2,qcm2);

% -------------------------------------------------------------
% Provide useful information about this patch.  Everything will be stored in
% the 'UserData' field of the patch handle p.
% --------------------------------------------------------------

% q data
userdata.q = qcm2;
userdata.qmin = min(min(qcm2));
userdata.qmax = max(max(qcm2));

% Spatial information
userdata.sdir = sdir;
userdata.sval = sval;

userdata.mappedgrid = mappedgrid;
userdata.manifold = manifold;

[vc{1:3}] = get_xyzlike(xc,yc,zc,sdir);
[ve{1:3}] = get_xyzlike(xe,ye,ze,sdir);
[vc_names{1:3}] = get_xyzlike('xc','yc','zc',sdir);
[ve_names{1:3}] = get_xyzlike('xe','ye','ze',sdir);
[mv_names{1:3}] = get_xyzlike('mx','my','mz',sdir);
[dv_names{1:3}] = get_xyzlike('dx','dy','dz',sdir);

for i = 1:3,
  userdata = setfield(userdata,(vc_names{i}),vc{i});
  userdata = setfield(userdata,(ve_names{i}), ve{i});
  userdata = setfield(userdata,(mv_names{i}),length(vc{i}));

  % If we are doing a manifold, then dz==1.
  dv = ve{i}(2) - ve{i}(1);
  userdata = setfield(userdata,(dv_names{i}),(dv == 0) + dv*(dv ~= 0));
end;

[vmin_names{1:3}] = get_xyzlike('xmin','ymin','zmin',sdir);
[vmax_names{1:3}] = get_xyzlike('xmax','ymax','zmax',sdir);
userdata = setfield(userdata,vmin_names{1},sval);
userdata = setfield(userdata,vmax_names{1},sval);
for i = 2:3,
  userdata = setfield(userdata,(vmin_names{i}), ve{i}(1));
  userdata = setfield(userdata,(vmax_names{i}), ve{i}(end));
end;

% Store Cartesian coodinates for use in mask_patches, and convert patch
% vertices to physical coordinates.
v = get(p,'Vertices');
userdata.cartCoords = v;
if (mappedgrid == 1 | manifold == 1)
  if (mappedgrid == 1)
    if (nargin('mapc2p') == 2)
      [v(:,1),v(:,2)] = mapc2p(v(:,1),v(:,2));
    else
      [v(:,1),v(:,2),v(:,3)] = mapc2p(v(:,1),v(:,2),v(:,3));
    end;
  end;
  if (manifold == 1)
    [v(:,1),v(:,2), v(:,3)] = mapc2m(v(:,1),v(:,2));
  end;
  set(p,'Vertices',v);
end;


% -------------------------------------------------------
% Now create and store some other graphic objects that are associated with
% this patch.
% -------------------------------------------------------
% Contour lines.
userdata.contourLines = [];
if (~isempty(contourlevels))
  c = contourc(yc_like,zc_like,qcm2,contourlevels,'k');
  userdata.contourLines = create_clines(c,sval,sdir,mappedgrid,manifold);
end;

% Mesh data for showing coarsened mesh later...
% userdata.mesh = create_mesh(sdir,sval,xe,ye,ze,mappedgrid,manifold);
% The mesh is created with showmesh.
userdata.mesh.xlines = [];
userdata.mesh.ylines = [];

% Lines at intersections of x,y,z planes - for 3d only.  These are actually
% created after all slices have been plotted.
userdata.grid_number = grid_number;  % For computing intersections
userdata.xyIntersect = [];
userdata.xzIntersect = [];
userdata.yzIntersect = [];

% Patch borders,
userdata.border  = create_border(sdir, sval, xe,ye,ze,mappedgrid,manifold);

% Set patch UserData.
set(p,'UserData',userdata);

% Set state
set(p,'Tag','on');
