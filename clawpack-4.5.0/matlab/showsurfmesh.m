function showsurfmesh(snum)

% SHOWSURFMESH shows mesh lines on isosurfaces.
%
%   SHOWSURFMESH(SNUM) shows mesh lines on surfaces specified in vector SNUM,
%   corresponding to entries in the ISOSURFVALUES parameter used to create
%   the isosurfaces.
%
%   See also HIDESURFMESH.

isurfaces = get_isosurfaces;

if (nargin == 0)
  snum = 1:length(isurfaces);
end;

for ns = 1:length(snum),
  n  = snum(ns);
  if (n < 0 | n > length(isurfaces))
    continue;
  end;
  isurfs = isurfaces{n};
  for l = 1:length(isurfs),
    isurf_vec = isurfs{l};
    for k = 1:length(isurf_vec),
      is = isurf_vec(k);
      set(is,'EdgeColor','k');
    end;
  end;
end;
