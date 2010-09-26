function setslicecolor(c,sdir,snum)

% SETSLICECOLOR sets the color of specified slices
%
%    SETSLICECOLOR(C,SDIR,SNUM) sets the color of slices corresponding to
%    xSliceCoords(SNUM), ySliceCoords(SNUM), or zSliceCoords(SNUM)
%    (depending on value of SDIR) to color C.  C should be a string variable
%    specifying a legitimate color value, or a 1x3 RGB vector.
%
%    SETSLICECOLOR(C,SDIR) sets all slices in direction SDIR to C.
%
%    SETSLICECOLOR(C) sets all slices to color C.
%
%    This is useful for creating black and white contour plots.
%
%    The color value may be either a legitimate character string, or an RGB
%    triple.
%
%    See also PLOT.



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
    for level = 1:length(slice),
      pvec = slice{level};
      set(pvec,'FaceColor',c);
    end;
  end;
end;
