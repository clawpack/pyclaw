function plotcubeedges = setPlotCubeEdges(plotcubeedges)

% setPlotCubeEdges shows or hides 3d amr cubes
%
%        PlotCubeEdges = setPlotCubeEdges(PlotCubeEdges) shows
%        amr cubes at level N if PlotCubeEdges(N) == 1 and hides the
%        cube otherwise.
%
%        By setting variable 'PlotCubeEdges' as the return argument, the effects
%        of hiding or showing the computational grid will carry over to the
%        next time Frame.
%
%        See also SHOWCUBES, HIDECUBES.

for level = 1:length(plotcubeedges)
  if (plotcubeedges(level) == 1)
    showcubes(level);
  else
    hidecubes(level);
  end;
end;
