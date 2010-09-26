function plotgridedges = setPlotGridEdges(plotgridedges)


% setPlotGridEdges shows or hides patch borders.
%
%        PlotGridEdges = setPlotGrid(PlotGridEdges) shows patch borders
%        of level N patches if PlotGridEdges(N) == 1 and hides the patch
%        borders at level N otherwise.
%
%        By setting variable 'PlotGridEdges' as the return argument,
%        the effects of hiding or showing the patch borders will
%        carry over to the next time Frame.
%
%        See also setPlotMesh, setPlotGrid, SHOWGRIDLINES, HIDEGRIDLINES.


for level = 1:length(plotgridedges)
  if (plotgridedges(level) == 1)
    showpatchborders(level);
  else
    hidepatchborders(level);
  end;
end;
