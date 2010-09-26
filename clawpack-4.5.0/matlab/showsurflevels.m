function showsurflevels(level)

% SHOWSURFLEVELS shows isosurfaces on specified levels.
%
%       SHOWSURFLEVELS(LEVEL) shows isosurface at amr levels specified in vector
% 	  LEVEL.
%
%       SHOWSURFLEVELS, by itself,  shows isosurfaces at all amr levels.
%
%       See also HIDESURFLEVELS.



isurfaces = get_isosurfaces;
for n = 1:length(isurfaces),
  isurfs = isurfaces{n};
  if (nargin == 0)
    level = 1:length(isurfs);
  end;
  for l = 1:length(level),
    isurf_vec = isurfs{level(l)};
    for k = 1:length(isurf_vec),
      is = isurf_vec(k);
      set(is,'Tag','on');
      set_isosurface_visibility(is,'on');
    end;
  end;
end;
