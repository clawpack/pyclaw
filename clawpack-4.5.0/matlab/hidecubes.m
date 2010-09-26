function hidecubes(level)

% HIDECUBES hides 3d patch borders
%
%   HIDECUBES(LEVEL) hides 3d amr cubes at levels specified by vector
%   LEVEL.
%
%   HIDECUBES, by itself, hides 3d amr cubes at all levels.
%
%   See SETPLOT3, SHOWCUBES, SHOWPATCHBORDERS, HIDEPATCHBORDERS.

cubes = get_cubes;

if (nargin == 0)
  level = 1:length(cubes);
end;

for l = 1:length(level),
  cube_vec = cubes{level(l)};
  for k = 1:length(cube_vec),
    cube = cube_vec(k);
    set(cube.c,'Visible','off');
  end;
end;
