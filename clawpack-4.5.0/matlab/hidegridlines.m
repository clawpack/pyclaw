function hidegridlines(level,sdir,snum)

% HIDEGRIDLINES hides computational grid lines
%
%     HIDEGRIDLINES(LEVEL) hides computational grid on levels specified
%     as in vector LEVEL.
%
%     Example :
%
%                hidegridlines(2:3);
%
%     hides grid lines on levels 2 and 3.
%
%     HIDEGRIDLINES, by itself, will hide computational grid on all levels.
%
%     See also SHOWGRIDLINES, SHOWMESH, HIDEMESH, SetPlotGrid.


if (nargin < 2)
  sdirs = {'x','y','z'};
else
  sdirs = {sdir};
end;

for idir = 1:length(sdirs),
  slices = get_slices(sdirs{idir});
  if (nargin < 3)
    snum = 1:length(slices);
  end;
  for ns = 1:length(snum),
    n = snum(ns);
    slice = slices{n};
    if (nargin == 0)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      set(pvec,'EdgeColor','none');
    end;
  end;
end;
