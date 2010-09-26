function set_isosurfaces(isurfaces)

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');

if (~strcmp(ftag,'AMRClawSlicePlot'))
  error('get_slices : Current figure does not contain slice data');
end;

amrplot = get(gcf,'UserData');

amrplot.isosurfaces = isurfaces;
set(gcf,'UserData',amrplot);
