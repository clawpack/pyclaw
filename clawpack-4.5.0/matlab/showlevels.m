function showlevels(level)

% SHOWLEVELS shows slice data
%
%       SHOWLEVELS(level) shows slice data at amr levels specified in vector
% 	  LEVEL.
%
%       SHOWLEVELS, by itself, shows slice data at all amr levels.
%
%       In 3d, SHOWLEVELS will show data on all slices, on specified
%       levels. For this reason, SHOWLEVELS is probably not very useful in
%       the 3d setting in its current implementation.
%
%       See also HIDELEVELS.


sdirs = {'x','y','z'};

for idir = 1:3,
  slices = get_slices(sdirs{idir});
  for n = 1:length(slices),
    slice = slices{n};
    if (nargin == 0)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      for k = 1:length(pvec),
	set(pvec(k),'Tag','on');
	% I don't set the visibility because I don't know which slices are
	% currently 'on'.
	set_patch_visibility(pvec(k),'on');
      end;
    end;
  end;
end;
