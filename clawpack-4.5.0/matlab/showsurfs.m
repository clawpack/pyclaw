function showsurfs(snum)

% SHOWSURFS shows isosurfaces.
%
%   SHOWSURFS(SNUM) shows isosurfaces in vector SNUM, where SNUM is a
%   vector of integers corresponding to entries in the IsosurfValues
%   parameter used to create the isosurfaces.
%
%   SHOWSURFS, by itself, shows all isosurfaces at all levels.
%
%   See also HIDESURFS, SURFLOOP.

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
  for level = 1:length(isurfs),
    isurf_vec = isurfs{level};
    for k = 1:length(isurf_vec),
      is = isurf_vec(k);
      set_isosurface_visibility(is,'on');
    end;
  end;
end;
