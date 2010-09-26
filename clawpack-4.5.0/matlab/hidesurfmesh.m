function hidesurfmesh(snum)

% HIDESURFMESH hides mesh lines on isosurfaces.
%
%   HIDESURFMESH(SNUM) hides mesh lines on surfaces specified in vector SNUM,
%   corresponding to entries in the ISOSURFVALUES parameter used to create
%   the isosurfaces.
%
%   See also SHOWSURFMESH.

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
      set(is,'EdgeColor','none');
    end;
  end;
end;
