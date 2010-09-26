function slices = get_slices(sdir)

% Internal matlab routine for Clawpack graphics.

ftag = get(gcf,'Tag');

if (~strcmp(ftag,'AMRClawSlicePlot'))
  error('get_slices : Current figure does not contain slice data');
end;

amrplot = get(gcf,'UserData');

slice_handles = amrplot.slices;

idir = findstr(lower(sdir),'xyz');
if (isempty(idir))
  error('get_slices : sdir must be equal to ''x'', ''y'', or ''z''');
end;

slices = slice_handles{idir};
