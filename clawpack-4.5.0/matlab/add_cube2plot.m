function add_cube2plot(xe,ye,ze,level,mappedgrid)

% Internal matlab routine for Clawpack graphics

% cube is a struct, whose field 'c' is a vector of handles to line objects.
cube = create_cube(xe,ye,ze,mappedgrid);

cubes = get_cubes;

cube_vec = cubes{level};
if (isempty(cube_vec))
  % the (end+1) trick doesn't work when trying to assign a struct to [].
  cube_vec = cube;
else
  cube_vec(end+1) = cube;
end;

cubes{level} = cube_vec;
set_cubes(cubes);
