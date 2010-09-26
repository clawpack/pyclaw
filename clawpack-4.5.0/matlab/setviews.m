% SETVIEWS sets view coordinates for pre-defined 3d viewing angles.
%
%     SETVIEWS has no arguments.
%
%     Once the setviews command has put viewing variables into the global
%     workspace, the following arguments are available to the view command:
%
%     vtop       - view plot from the top
%     vbot       - view plot from the bottom
%     vleft      - view plot from the left
%     vright     - view plot from the right
%     vback      - view plot from the back
%     vfront     - view plot from the front.
%
%     xslice     - view an x-slice. Equivalent to vright.
%     yslice     - view a y-slice. Equivalent to vfront.
%     zslice     - view a z-slice. Equivalent to vtop.
%
%     Example :
%                 setviews;
%                 surf(peaks);
%                 view(vtop);  % view 'peaks' plot from the top.
%
%     See also VIEW.

global vtop  vbot vfront  vback  vleft  vright  xslice  yslice zslice;

vtop = [0,90];
vbot = [0,-90];
vfront = [0,0];
vback = [180,0];
vleft = [-90,0];
vright = [90,0];


xslice = vright;
yslice = vfront;
zslice = vtop;
