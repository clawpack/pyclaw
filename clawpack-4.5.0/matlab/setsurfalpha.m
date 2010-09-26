function setsurfalpha(a,snum)

% SETSURFALPHA sets the alpha (transparency) value for isosurfaces
%
%   SETSURFALPHA(A,SNUM) sets the alpha value to A for isosurfaces
%   specified in vector SNUM, where the entries of SNUM are integers
%   specifying isosurface number corresponding to entries in ISOSURFVALUES.
%   A should be a value between 0 and 1.
%
%   SETSURFALPHA(A) sets all isosurfaces to alpha value A.
%
%   SETSURFALPHA, by itself, sets the alpha value to 1.
%
%   SETSURFALPHA sets the 'FaceAlpha' property of the isosurface patch.
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
%   See also SETSURFCOLOR, SETOPENGL, OPENGL, PATCH, ISOSURFACE.


rstr = get(gcf,'Renderer');
if (strcmp(rstr,'OpenGL') == 0)
  return;
end;

a = min([max([a,0]),1]);

isurfaces = get_isosurfaces;

if (nargin < 2)
  snum = 1:length(isurfaces);
  if (nargin < 1)
    a = 1;
  end;
end;

for ns = 1:length(snum),
  n = snum(ns);
  if (n < 1 | n > length(isurfaces))
    continue;
  end;
  isurfs = isurfaces{n};
  for level = 1:length(isurfs),
    isurf_vec = isurfs{level};
    for k = 1:length(isurf_vec),
      is = isurf_vec(k);
      set(is,'FaceAlpha',a);
      udata = get(is,'UserData');
      udata.alpha = a;
      set(is,'UserData',udata);
    end;
  end;
end;
