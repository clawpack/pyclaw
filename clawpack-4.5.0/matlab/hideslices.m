function hideslices(sdir,snum)

% HIDESLICES hides 3d slices
%
%   HIDESLICES(SDIR) hides all slices in direction SDIR.  SDIR can
%   be equal to 'x', 'y', or 'z'.
%
%   HIDESLICES, by itself, hides all slices in all directions.
%
%   See also SHOWSLICES, SLICELOOP.
%

if (nargin == 0)
  sdirs = {'x', 'y', 'z'};
else
  sdirs = {sdir};
end;

for idir = 1:length(sdirs),
  slices = get_slices(sdirs{idir});
  if (nargin <= 1)
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
      for k = 1:length(pvec),
	set_patch_visibility(pvec(k),'off');
      end; % Patches loop
    end; % level loop
  end;
end;
