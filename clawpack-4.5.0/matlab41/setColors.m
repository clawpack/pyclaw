function setColors(p,q)


% Generic way to set colors.
% The 'q' data will be used to determine scaled colors in the
% current colormap.  The user may specify their own version of this
% to get a different colormap.

set(p,'FaceVertexCData',q);
set(p,'FaceColor','flat');
