function h = projectBorder(xlow,xhigh,ylow,yhigh,zlow,zhigh,sdir,sval)

% This projects the border onto an underlying slice slayer.
% it is assumed that the cube described by xlow,xhigh, etc actually
% intersects the slice.

% Now plot boxes directly on slices.
if (strcmp(sdir,'x'))
  xlow = sval;
  xhigh = sval;
elseif (strcmp(sdir,'y'))
  ylow = sval;
  yhigh = sval;
elseif (strcmp(sdir,'z'))
  zlow = sval;
  zhigh = sval;
end;

h = plotCube(xlow,xhigh,ylow,yhigh,zlow,zhigh);
