function isurf = get_isosurfaces()

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');

if (~strcmp(ftag,'AMRClawSlicePlot'))
  error('get_isosurfaces : Current figure does not contain slice data');
end;

amrplot = get(gcf,'UserData');

isurf = amrplot.isosurfaces;
