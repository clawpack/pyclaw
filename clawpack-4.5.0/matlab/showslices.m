function showslices(sdir,snum)

% SHOWSLICES makes 3d slices visible
%
%   SHOWSLICES(SDIR) shows all slices in direction SDIR.  SDIR can
%   be equal to 'x', 'y', or 'z'.
%
%   SHOWSLICES(SDIR,SNUM) shows slices in direction SDIR specified by slices
%   numbers given in vector SNUM.
%
%   SHOWSLICES, by itself, shows all slices in all directions.
%
%   In order for a slice to be visible, it must have been created by setting
%   xSliceCoords, ySliceCoords or zSliceCoords in setplot3.m
%
%   See also HIDESLICES, SLICELOOP, SETPLOT.

if (nargin < 1)
  sdirs = {'x', 'y', 'z'};
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
      for k = 1:length(pvec),
	set_patch_visibility(pvec(k),'on');
      end; % Patches loop
    end; % level loop
  end;
end;
