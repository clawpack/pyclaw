function hidecontourlines(sdir,snum)

% HIDECONTOURLINES hides contourlines.
%
%     HIDECONTOURLINES has no arguments.
%
%     HIDECONTOURLINES hides all contours specified in user defined parameter
%     ContourValues, set in SETPLOT2 or SETPLOT3.  All contour lines, on all
%     levels and slices, in all directions will be hidden.
%
%     If ContourValues is initially empty, or doesn't exist, no contours lines
%     will have been created, and HIDECONTOURLINES has no effect.
%
%     See also SHOWCONTOURLINES, SETPLOT2, SETPLOT3.
%

if (nargin < 1)
  sdirs = {'x', 'y','z'};
else
  sdirs = {sdir};
end

for idir = 1:length(sdirs),
  slices = get_slices(sdirs{idir});
  if (nargin < 2)
    snum = 1:length(slices);
  end;
  for ns = 1:length(snum),
    n = snum(ns);
    slice = slices{n};
    for level = 1:length(slice),
      pvec = slice{level};
      for k = 1:length(pvec),
	p = pvec(k);
	udata = get(p,'UserData');
	set(udata.contourLines,'Tag','off');
	set_cline_visibility(pvec(k));
      end;
    end;
  end;
end;
