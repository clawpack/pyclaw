function showgridlines(level,sdir,snum)

% SHOWGRIDLINES shows computational grid lines
%
%     SHOWGRIDLINES(LEVEL,SDIR,SNUM) shows gridlines at levels in vector
%     LEVEL, in a direction specified by SDIR, ('x','y' or 'z') on slices
%     specified in vector SNUM, where the entries of SNUM are integers
%     corresponding to slices in user specified parameters xSliceCoorrds,
%     ySliceCoords or zSliceCoords.
%
%     SHOWGRIDLINES(LEVEL) shows computational grid on levels specifed in
%     vector LEVEL.  Gridlines in all directions, on all slices will be shown.
%
%     SHOWGRIDLINES, by itself, will show computational grid on all levels.
%
%     Example :
%
%                showgridlines([1 2]);
%
%     shows grid lines on levels 1 and 2.
%
%                showgridlines(1,'x',5)
%
%     shows gridlines at level 1 on slice corresponding to xSliceCoords(5).
%
%     See also HIDEGRIDLINES, SHOWMESH, HIDEMESH, SetPlotGrid.


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
    if (n < 1 | n > length(slices))
      continue;
    end;
    slice = slices{n};
    if (nargin == 0)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      set(pvec,'EdgeColor','k');
    end;
  end;
end;
