function clear_amrplot()

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');
if (strcmp(ftag,'AMRClawSlicePlot') == 0)
  % nothing to clear in current window
  return;
end;

sdir = {'x', 'y', 'z'};

for idir = 1:3,
  slices = get_slices(sdir{idir});
  for n = 1:length(slices),
    slice = slices{n};
    for level = 1:length(slice),
      pvec = slice{level};
      for k = 1:length(pvec), % Loop over patches at level 'level'
	if (ishandle(pvec(k)))
	  delete_patch(pvec(k));
	end;
      end;
    end;
  end;
end;

cubes = get_cubes;
for level = 1:length(cubes),
  cube_vec = cubes{level};
  for k = 1:length(cube_vec),
    cube = cube_vec(k);
    if (ishandle(cube.c))
      delete(cube.c);
    end;
  end;
end;

isurfaces = get_isosurfaces;
for n = 1:length(isurfaces),
  isurf = isurfaces{n};
  for level = 1:length(isurf),
    isurf_vec = isurf{level};
    for k = 1:length(isurf_vec),
      is = isurf_vec(k);
      if (ishandle(is))
	delete(is);
      end;
    end;
  end;
end;

lineplot = get_lines;
for level = 1:length(lineplot),
  lvec = lineplot{level};
  for k = 1:length(lvec),
    lp = lvec(k);
    if (ishandle(lp))
      delete(lp);
    end;
  end;
end;


set(gcf,'Tag','');
set(gcf,'UserData',[]);
