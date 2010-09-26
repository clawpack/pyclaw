function hidepatchborders(level)

% HIDEPATCHBORDERS hides patch borders
%
%    HIDEPATCHBORDERS(level) hides patch borders for all patches at levels
%    specified in vector LEVEL.
%
%    HIDEPATCHBORDERS, by itself,  hides patch borders on all amr
%    patches for all levels
%
%    See also SHOWPATCHBORDERS, SETPLOTGRIDEDGES


sdir = {'x','y','z'};

for idir = 1:3,
  slices = get_slices(sdir{idir});
  for n = 1:length(slices),
    slice = slices{n};
    if (nargin == 0)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      for k = 1:length(pvec),
	p = pvec(k);
	udata = get(p,'UserData');
	set(udata.border,'Tag','off');
	set_pborder_visibility(p);
      end;
    end;
  end;
end;
