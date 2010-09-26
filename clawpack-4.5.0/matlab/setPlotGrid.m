function plotgrid = setPlotGrid(plotgrid)

% setPlotGrid shows or hides computational grid.
%
%        PlotGrid = setPlotGrid(PlotGrid) shows computational grid at level
%        N if PlotGrid(N) == 1 and hides the grid otherwise.
%
%        By setting variable 'PlotGrid' as the return argument, the effects
%        of hiding or showing the computational grid will carry over to the
%        next time Frame.
%
%        See also setPlotGridEdges, setPlotData, SHOWGRIDLINES,
%        HIDEGRIDLINES.


for level = 1:length(plotgrid)
  if (plotgrid(level) == 1)
    showgridlines(level);
  else
    hidegridlines(level);
  end;
end;
