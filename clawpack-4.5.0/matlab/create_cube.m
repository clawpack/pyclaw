function cube = create_cube(xe,ye,ze,mappedgrid);

% Internal Matlab routine for Clawpack graphics.

manifold = 0;

sdirs = {'x', 'y', 'z'};
ve = {xe, ye, ze};

bc = cell(6,1);
for idir = 1:3,
  sval = ve{idir}(1);
  bc{2*idir - 1} = create_border(sdirs{idir},sval,xe,ye,ze,mappedgrid,manifold);

  sval = ve{idir}(end);
  bc{2*idir} = create_border(sdirs{idir},sval,xe,ye,ze,mappedgrid,manifold);
end;

c = [bc{:}];

set(c,'Visible','off');

cube.c = c;  % c is a vector of handles to line objects.
