function hidesurfs(snum)

% HIDESURFS hides isosurfaces
%
%   HIDESURFS(SNUM) hides isosurfaces in vector SNUM, where SNUM is a
%   vector of integers corresponding to entries in the IsosurfValues
%   parameter used to create the isosurfaces.
%
%   HIDESURFS, by itself, hides all isosurfaces at all levels.
%
%   See SHOWSURFS, SURFLOOP, SETPLOT.

isurfaces = get_isosurfaces;

if (nargin == 0)
  snum = 1:length(isurfaces);
end;

for ns = 1:length(snum),
  n = snum(ns);
  if (n < 1 | n > length(isurfaces))
    continue;
  end;
  isurfs = isurfaces{n};
  for l = 1:length(isurfs),
    isurf_vec = isurfs{l};
    for k = 1:length(isurf_vec),
      is = isurf_vec(k);
      set_isosurface_visibility(is,'off');
    end;
  end;
end;
