function showcontourlines(sdir,snum)

% SHOWCONTOURLINES shows contourlines.
%
%     SHOWCONTOURLINES(SDIR,SNUM) shows contour lines on slices
%     xSliceCoords(SNUM), ySliceCoords(SNUM) or zSliceCoords(SNUM),
%     depending on value SDIR (='x','y','z').
%
%     SHOWCONTOURLINES(SDIR) shows contour lines on all slices in direction
%     SNUM.
%
%     SHOWCONTOURLINES, by itself, shows contour lines on all slices in all
%     directions.
%
%     SHOWCONTOURLINES shows contour lines that have been specified in the
%     ContourValues plotting parameter in SETPLOT2 or SETPLOT3.
%
%     If ContourValues is initially empty, or doesn't exist, no contours lines
%     will have been created, and SHOWCONTOURLINES has no effect.
%
%     See also HIDECONTOURLINES, SETPLOT2, SETPLOT3.

if (nargin < 1)
  sdirs = {'x', 'y','z'};
else
  sdirs = {sdir};
end;
for idir = 1:length(sdirs),
  slices = get_slices(sdirs{idir});
  if (nargin < 2)
    snum = 1:length(slices);
  end;
  for ns = 1:length(snum),
    n = snum(ns);
    if (n < 1 | n > length(slices))
      continue;
    end;
    slice = slices{n};
    for level = 1:length(slice),
      pvec = slice{level};
      for k = 1:length(pvec);
	p = pvec(k);
	udata = get(p,'UserData');
	set(udata.contourLines,'Tag','on');
	set_cline_visibility(pvec(k));
      end;
    end;
  end;
end;
