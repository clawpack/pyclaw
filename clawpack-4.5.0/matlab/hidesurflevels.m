function hidesurflevels(level)

% HIDESURFLEVELS hides isosurfaces on specified levels.
%
%       HIDESURFLEVELS(level) hides isosurface at amr levels specified in vector
% 	  LEVEL.
%
%       HIDESURFLEVELS, by itself,  hides isosurfaces data at all amr levels.
%
%       See also SHOWSURFLEVELS.



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
      set(is,'Tag','off');
      set_isosurface_visibility(is,'off');
    end;
  end;
end;
