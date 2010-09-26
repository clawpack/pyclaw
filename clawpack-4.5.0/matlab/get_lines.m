function lh = get_lines()

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');

if (~strcmp(ftag,'AMRClawSlicePlot'))
  error('get_lines : Current figure does not contain slice data');
end;

amrplot = get(gcf,'UserData');

lh = amrplot.lines;
