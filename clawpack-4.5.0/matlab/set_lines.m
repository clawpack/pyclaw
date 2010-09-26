function set_lines(amrlines);

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');

if (~strcmp(ftag,'AMRClawSlicePlot'))
  error('get_slices : Current figure does not contain slice data');
end;

amrplot = get(gcf,'UserData');

amrplot.lines = amrlines;
set(gcf,'UserData',amrplot);
