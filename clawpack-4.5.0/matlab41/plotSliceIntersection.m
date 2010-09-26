function plotSliceIntersection(ax,bx,ay,by,az,bz,sdir1,sdir2,sval1,sval2);

% This plots the line represented by the intersection of the two slices, s1
% and s2.

px = [];
py = [];
pz = [];


if (strcmp(sdir1,'x'))
  px = [sval1 sval1];
elseif (strcmp(sdir1,'y'))
  py = [sval1 sval1];
elseif (strcmp(sdir1,'z'))
  pz = [sval1 sval1];
end;

if (strcmp(sdir2,'x'))
  px = [sval2 sval2];
elseif (strcmp(sdir2,'y'))
  py = [sval2 sval2];
elseif (strcmp(sdir2,'z'))
  pz = [sval2 sval2];
end;


if (isempty(px))
  px = [ax bx];
elseif (isempty(py))
  py = [ay by];
elseif (isempty(pz))
  pz = [az bz];
end;

h = line('XData',px,'YData',py,'ZData',pz,'Color','k');
