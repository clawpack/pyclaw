function setsurfcolor(c,snum)

% SETSURFCOLOR sets the color for a specified isosurface
%
%   SETSURFCOLOR(A,SNUM) sets the color value to C for isosurfaces
%   specified in vector SNUM, where the entries of SNUM are integers
%   specifying isosurface number corresponding to entries in ISOSURFVALUES.
%   C should be a legitimate color string ('b','r',etc) or a RGB triple.
%
%   SETSURFCOLOR(C) sets all isosurfaces to color C.
%
%   SETSURFCOLOR, by itself, sets the isosurfaces to 'b'.
%
%   SETSURFCOLOR sets the 'FaceCOLOR' property of the isosurface patch.
%
%   See also SETSURFALPHA, PATCH, ISOSURFACE.


isurfaces = get_isosurfaces;

if (nargin < 2)
  snum = 1:length(isurfaces);
  if (nargin < 1)
    c = 'b';
  end;
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
      set(is,'FaceColor',c);
    end;
  end;
end;
