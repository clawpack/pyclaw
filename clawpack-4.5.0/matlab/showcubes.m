function showcubes(level,c,lw)

% SHOWCUBES shows 3d patch borders
%
%   SHOWCUBES(LEVEL,C,LW) shows cubes at levels specified by vector
%   LEVEL, and sets the edge color to C, and linewidth to LW.
%
%   SHOWCUBES(LEVEL,C) uses default linewidth
%
%   SHOWCUBES(LEVEL) shows 3d amr cubes at levels specified by vector
%   LEVEL.  Color is 'k' and uses default linewidth.
%
%   SHOWCUBES, by itself, shows 3d amr cubes at all levels.
%
%   See SETPLOT3, HIDECUBES, SHOWPATCHBORDERS, HIDEPATCHBORDERS.

cubes = get_cubes;

if (nargin < 3)
  lw = 1;
  if (nargin < 2)
    c = 'k';
    if nargin < 1
      level = 1:length(cubes);
    end;
  end;
end;

for l = 1:length(level),
  cube_vec = cubes{level(l)};
  for k = 1:length(cube_vec),
    cube = cube_vec(k);
    set(cube.c,'Visible','on','Color',c,'LineWidth',lw);
  end;
end;
