function create_amrplot(MaxLevels,xscoords, yscoords, ...
    zscoords,isosurfvalues)

% Internal matlab routine for Clawpack graphics.

if (nargin < 5)
  if (nargin <= 4)
    isosurfvalues = [];
    if (nargin == 1)
      xscoords = [];
      yscoords = [];
      zscoords = [];
    end;
  end;
end;

sliceCoords = {xscoords, yscoords, zscoords};
sdir = {'x', 'y', 'z'};

for idir = 1:3,
  n = length(sliceCoords{idir}); % Number of slices in this direction
  if (n > 0)
    % Fill up each slice with MaxLevel empty arrays;  These arrays will then
    % be filled with patches at each level.
    [slice_handles{idir}{1:n,1}] = deal(cell(1,MaxLevels));
  else
    % There are no slices in this direction.
    slice_handles{idir} = {};
  end;
end;

% Set various things in the amr plot.
amrplot.slices = slice_handles;

% Cubes
[amrplot.cubes{1:MaxLevels}] = deal([]);

% Isosurfaces
n = length(isosurfvalues);
if (n > 0)
  [amrplot.isosurfaces{1:n}] = deal(cell(1,MaxLevels));
else
  amrplot.isosurfaces = {};
end;

% Scatter plot
amrplot.lines = cell(MaxLevels,1);

set(gcf,'UserData',amrplot);
set(gcf,'Tag','AMRClawSlicePlot');
