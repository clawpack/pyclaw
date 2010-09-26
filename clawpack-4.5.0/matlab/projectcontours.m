function p_out = projectcontours(sdir,sval)
%
% PROJECTCONTOURS projects contours lines from a manifold to a specifed x,
%                 y, or z plane.
%
%     PROJECTCONTOURS(SDIR,SVAL) projects contours from manifold to plane at
%     constant SDIR at value SVAL.   The axes are automatically adjusted
%     in the direction SDIR so that the newly created contour lines are
%     visible.
%
%     P = PROJECTCONTOURS(SDIR,SVAL) returns a patch object whose dimensions
%     are specified by the current axis limits.  Initially the 'FaceColor'
%     property of this patch object is set to 'w', and the 'EdgeColor' is
%     set to 'k'.  The user may change these properties with the following
%     commands:
%
%     Example :
%           p = projectcontours('z',-1);  % projects contours from manifold
%                                         %    plane at z = -1.
%           set(p,'EdgeColor','none');    % Edge of patch containing new
%                                         %    contours is made invisible.
%           set(p,'FaceColor','r');       % Set patch face color to red.
%
%
%     NOTE : If no contours were created on original manifold, this command
%     does nothing.  To create contours, set ContourValues to scalar (number
%     of contour lines) or vector of contour values.
%
%     This command works only with 2d plots.
%
%     See also PATCH, MANIFOLD, SHOWCONTOURLINES, HIDECONTOURLINES.


slices = get_slices('z');
if (length(slices) ~= 1)
  error('projectcontours : There should be exactly one slice in the z direction');
end;

slice = slices{1};  % assume one slice in z direction.

xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
zlim = get(gca,'ZLim');

[xlim_like, ylim_like, zlim_like] = get_xyzlike(xlim,ylim,zlim,sdir);
[yem2_like,zem2_like] = meshgrid(ylim_like, zlim_like);
xem2_like = 0*yem2_like + sval;
[xem2, yem2, zem2] = get_xyz(xem2_like, yem2_like, zem2_like,sdir);

% Create new patch and set color to white.
p_new = patch(surf2patch(xem2,yem2, zem2));

set(p_new,'FaceColor','w');
set(p_new,'EdgeColor','k');

clines_new = [];

for level = 1:length(slice), % loop over levels
  pvec = slice{level};  % patches at level 'level'
  for k = 1:length(pvec),
    p = pvec(k);
    udata = get(p,'UserData');

    % now project the contourlines.
    clines = udata.contourLines;
    hdl = [];
    for i = 1:length(clines),
      xdata = get(clines(i),'XData');
      ydata = get(clines(i),'YData');
      zdata = get(clines(i),'ZData');
      [xdata_like, ydata_like, zdata_like] = get_xyzlike(xdata,ydata,zdata,sdir);
      xdata_like = 0*xdata_like + sval;
      [xdata,ydata,zdata] = get_xyz(xdata_like, ydata_like, zdata_like,sdir);
      hdl(i) = line('XData',xdata,'YData',ydata,'ZData',zdata);
      udata_cline = get(clines(i),'UserData');
      udata_h.cval = udata_cline.cval;
      set(hdl(i),'UserData',udata_h);
    end;
    clines_new = [clines_new hdl];
  end;
end;

udata_new.contourLines = clines_new;

% Set user data for new patch.
set(p_new,'UserData',udata_new);

% Change limits of axis so that we can see the projected contour lines.
prop = sprintf('%sLimMode',upper(sdir));
set(gca,prop,'auto');

if (nargout == 1)
  p_out = p_new;
end;
