function reducesurf(r,snum,level)

% REDUCESURF reduces the number of faces in specified isosurfaces
%
%   REDUCESURF(R,SNUM,LEVEL) reduces isosurfaces specified by entries in
%   vector SNUM, corresponding to values set in parameter ISOSURFVALUES (set
%   in setplot3.m) at levels specified in vector LEVEL.  Isosurface is
%   reduced from original isosurface by factor R, where 0 < R <= 1.
%
%   REDUCESURF(R,SNUM) reduces those isosurfaces specifed by entries in SNUM
%   at all levels.
%
%   REDUCESURF(R) reduces all isosurfaces at all levels by factor
%   R.
%
%   REDUCESURF, by itself, reduces all isosurfaces, at all levels by factor
%   R = 0.5.
%
%   Reducing the isosurface is useful for visualizing the isosurface mesh,
%   which, if shown at full resolution, maybe too dense to be of any use.
%
%   Use REDUCESURF(1) to restore original isosurface computed from data.
%
%   REDUCESURF uses the Matlab function REDUCEPATCH to reduce the isosurface.
%
%   See also SHOWSURFMESH, HIDESURFMESH, REDUCEPATCH.

isurfaces = get_isosurfaces;

if (nargin <= 1)
  if (nargin == 0)
    r = 0.5;
  end;
  snum = 1:length(isurfaces);
end;

lightprops = {'FaceLighting',...
	     'EdgeLighting',...
	     'BackFaceLighting',...
	     'AmbientStrength',...
	     'DiffuseStrength',...
	     'SpecularStrength',...
	     'SpecularExponent',...
	     'SpecularColorReflectance'};

for ns = 1:length(snum),
  n = snum(ns);
  if (n < 1 | n > length(isurfaces))
    continue;
  end;
  isurfs = isurfaces{n};
  if (nargin <= 2)
    level = 1:length(isurfs);
  end;
  for l = 1:length(level),
    isurf_vec = isurfs{level(l)};
    for k = 1:length(isurf_vec),
      is = isurf_vec(k);
      udata = get(is,'UserData');
      ivalue = udata.value;
      icolor = udata.color;
      ialpha = udata.alpha;
      mappedgrid = udata.mappedgrid;

      % Create new isosurface and reduce by factor r.
      is_new = create_isosurface(udata.xc,udata.yc,udata.zc,udata.q,...
	  r,ivalue,icolor,ialpha,mappedgrid);

      % Restore edgecolor of new patch.
      set(is_new,'EdgeColor',get(is,'EdgeColor'));

      % Restore edgecolor of new patch.
      set(is_new,'FaceColor',get(is,'FaceColor'));

      % Restore current lighting properties to new isosurface.
      for i = 1:length(lightprops),
	lp = get(is,lightprops{i});
	set(is_new,lightprops{i},lp);
      end;

      % Restore visiblity of this isosurface patch
      set(is_new,'Visible',get(is,'Visible'));

      isurfaces{n}{level(l)}(k) = is_new;
      delete(is);
    end;
  end;
end;

set_isosurfaces(isurfaces);
