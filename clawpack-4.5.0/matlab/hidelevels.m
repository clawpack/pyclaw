function hidelevels(level)

% HIDELEVELS hides slice data
%
%       HIDELEVELS, by itself,  hides slice data at all amr levels.
%
%       HIDELEVELS(level) hides slice data at amr levels specified in vector
% 	  LEVEL.
%
%       Note that HIDELEVELS is probably only useful for 2d plots for which
%       Manifold == 0.  If Manifold == 1, or 3d slices are being plotted,
%       patches on level k mask out regions on level k-1.  While HIDELEVEL(k)
%       will hide the patches on level k, data on level k-1 is not shown
%       in the regions occupied by level k patches.
%
%       See also SHOWLEVELS.



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
	set(pvec(k),'Tag','off');
	set_patch_visibility(pvec(k),'off');
      end;
    end;
  end;
end;
