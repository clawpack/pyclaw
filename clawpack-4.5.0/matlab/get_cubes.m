function c = get_cubes()

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');

if (~strcmp(ftag,'AMRClawSlicePlot'))
  error('get_cubes : Current figure does not contain slice data');
end;

amrplot = get(gcf,'UserData');

c = amrplot.cubes;
