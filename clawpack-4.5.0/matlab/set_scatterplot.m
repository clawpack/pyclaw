function set_scatterplot(splot)

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');

if (~strcmp(ftag,'AMRClawSlicePlot'))
  error('get_slices : Current figure does not contain slice data');
end;

amrplot = get(gcf,'UserData');

amrplot.scatterplot = splot;
set(gcf,'UserData',amrplot);
