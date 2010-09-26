function setopengl()

% SETOPENGL sets the graphics renderer to OpenGL.
%
%    SETOPENGL sets the current graphics renderer to OpenGL, if this
%    renderer is available on your system, and returns a warning otherwise.
%
%    Many 3d graphics features are rendered much better using OpenGL, so
%    it is suggested that it be used whenever possible.
%
%    SETOPENGL sets the 'Renderer' property of the current figure to
%    'OpenGL'.
%
%    Some graphics hardware may not support the OpenGL renderer, or may
%    cause certain Matlab routines to crash. For this reason, the choice
%    of whether to use OpenGL is left to the user.
%
%    See also OPENGL, FIGURE, GCF.

rset = set(gcf,'Renderer');
found_opengl = 0;
for i = 1:length(rset),
  if (strcmp(rset(i),'OpenGL'))
    found_opengl = 1;
    break;
  end;
end;

if (~found_opengl)
  disp('*** Warning : The OpenGL renderer is not available on your system.');
else
  set(gcf,'Renderer','OpenGL');
end;
