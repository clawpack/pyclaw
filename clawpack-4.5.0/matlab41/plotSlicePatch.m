function [p,h] = plotSlicePatch(xc,yc,zc,xe,ye,ze,wcm,sdir,sval,plotedges,contourLevels)

% Note that it is assumed that the cube, described by xe,ye,ze intersects
% the slice.

if (strcmp(sdir,'x'))
  % Get x slice.
  [yem,zem] = meshgrid(ye,ze);
  xem = 0*yem + sval;
elseif (strcmp(sdir,'y'))
  [xem,zem] = meshgrid(xe,ze);
  yem = 0*xem + sval;
elseif (strcmp(sdir,'z'))
  [xem,yem] = meshgrid(xe,ye);
  zem = 0*xem + sval;
end;

% Get patch vertices and faces from surface data.
fv = surf2patch(xem,yem,zem);

% Get centers of faces so that we can interpolate values from wcm to the
% slice using the location of the centers.

faces = fv.faces;
vertices = fv.vertices;

vx = vertices(:,1);
vy = vertices(:,2);
vz = vertices(:,3);
fc(:,1) = sum(vx(faces'))'/4;
fc(:,2) = sum(vy(faces'))'/4;
fc(:,3) = sum(vz(faces'))'/4;

% Now clamp face locations to the endpoints of xc,yc,zc. Otherwise,
% interp3 will return NaN's in locations where the data spills over the
% limits defined by xcm,ycm,zcm

% NOTE : If a slice is within dx/2 (or dy/2 or dz/2) of a patch boundary, we
% are actually going to use data at the center point as our 'extrapolated'
% data. This is probably okay, since at least visually, we'll get something
% that looks good.  We could possible detect this situation at the
% computational boundaries, and insist that slice values be within the range
% of the centered data, but this probably isn't necessary.  For interior
% patches, we could interpolate from coarser grids, but this would be a
% headache.  So we'll just live with this for now... In any case, the data
% will always appear to be plotted at the specified slice value,even if
% technically, that slice value falls outside of the range of

fc(find(fc(:,1) < xc(1)),1) = xc(1);
fc(find(fc(:,1) > xc(end)),1) = xc(end);

fc(find(fc(:,2) < yc(1)),2) = yc(1);
fc(find(fc(:,2) > yc(end)),2) = yc(end);

fc(find(fc(:,3) < zc(1)),3) = zc(1);
fc(find(fc(:,3) > zc(end)),3) = zc(end);


% Interpolate values to the centers of the faces
[xcm,ycm,zcm] = meshgrid(xc,yc,zc);

qslice = interp3(xcm,ycm,zcm,wcm,fc(:,1),fc(:,2),fc(:,3),'*linear');

p = patch(fv);

setColors(p,qslice);

user_data.sdir = sdir;
user_data.sval = sval;
user_data.qmin = min(qslice);
user_data.qmax = max(qslice);

if (strcmp(sdir,'x'))
  user_data.xmax = sval;
  user_data.xmin = sval;
  user_data.ymax = max(ye);
  user_data.ymin = min(ye);
  user_data.zmax = max(ze);
  user_data.zmin = min(ze);
elseif (strcmp(sdir,'y'))
  user_data.xmax = max(xe);
  user_data.xmin = min(xe);
  user_data.ymax = sval;
  user_data.ymin = sval;
  user_data.zmax = max(ze);
  user_data.zmin = min(ze);
elseif (strcmp(sdir,'z'))
  user_data.xmax = max(xe);
  user_data.xmin = min(xe);
  user_data.ymax = max(ye);
  user_data.ymin = min(ye);
  user_data.zmax = sval;
  user_data.zmin = sval;
end;

set(p,'UserData',user_data);

if (plotedges == 0)
  set(p,'EdgeColor','none');
end;

h = 0;
if (~isempty(contourLevels))
  h = plotSliceContour(xc,yc,zc,sval,sdir,xcm,ycm,zcm,wcm,contourLevels);
end;
