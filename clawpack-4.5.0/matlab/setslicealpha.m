function setslicealpha(a,sdir,snum)

% SETSLICEALPHA sets the alpha (transparency) value of specified slices
%
%    SETSLICECOLOR(A,SDIR,SNUM) sets the alpha value of slices corresponding to
%    xSliceCoords(SNUM), ySliceCoords(SNUM), or zSliceCoords(SNUM)
%    (depending on value of SDIR) to alpha value A.  A should be a value
%    between 0 and 1.  SDIR should be set to 'x', 'y', 'z'.
%
%    SETSLICEALPHA(A,SDIR) sets alpha value of all slices in direction
%    SDIR.
%
%    SETSLICEALPHA(A) sets all slices to alpha value A.
%
%    Transparency will only be set if the OpenGL renderer is the current
%    figure renderer.  To see if you have this renderer available on your
%    system, use the command
%
%           set(gcf,'Renderer')
%
%    to get list of available renderers for your system. If you have the OpenGL
%    renderer, you can set it using
%
%           set(gcf,'Renderer','OpenGL');
%
%    Or, simply use the command SETOPENGL.  This will set the renderer to
%    OpenGL if you have it on your system, and report a warning otherwise.
%
%    If the OpenGL Renderer is not set, this command has no effect.
%
%    See also OPENGL, SETOPENGL.

rstr = get(gcf,'Renderer');
if (strcmp(rstr,'OpenGL') == 0)
  return;
end;

if (nargin < 2)
  sdirs = {'x','y','z'};
  if (nargin < 1)
    a = 1.0;
  end;
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
    slice = slices{n};
    for level = 1:length(slice),
      pvec = slice{level};
      set(pvec,'FaceAlpha',a);
    end;
  end;
end;
