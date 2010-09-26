function showpatchborders(level,sdir,snum)

% SHOWPATCHBORDERS shows amr patch borders
%
%    SHOWPATCHBORDERS(level) shows patch borders for all patches at levels
%    specified by vector LEVEL.
%
%    SHOWPATCHBORDERS, by itself shows patch borders on all amr patches
%    for all levels
%
%    See also HIDEPATCHBORDERS, SETPLOTGRIDEDGES.
%

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
    if (nargin < 1)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      for k = 1:length(pvec),
	p = pvec(k);
	udata = get(p,'UserData');
	set(udata.border,'Tag','on');
	set_pborder_visibility(p);
      end;
    end;
  end;
end;
