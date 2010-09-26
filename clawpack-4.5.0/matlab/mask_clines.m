function mask_clines(p,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh)

% Internal matlab routine for Clawpack graphics.

udata = get(p,'UserData');

s = 1e-5;

[xlow_like, ylow_like, zlow_like]    = get_xyzlike(xlow,ylow,zlow,sdir);
[xhigh_like, yhigh_like, zhigh_like] = get_xyzlike(xhigh,yhigh,zhigh,sdir);

% Now mask any contour lines
ch = udata.contourLines;
for j = 1:length(ch),  % loop over contour lines on this patch
  udata = get(ch(j),'UserData');
  vc = udata.cartCoords;
  [xdata_like,ydata_like, zdata_like]  = get_xyzlike(vc(:,1),vc(:,2), vc(:,3),sdir);

  nan_mask = ydata_like > ylow_like+s & ydata_like < yhigh_like-s & ...
      zdata_like > zlow_like+s & zdata_like < zhigh_like-s;

  % Mask out values in physical data
  xdata = get(ch(j),'XData');
  ydata = get(ch(j),'YData');
  zdata = get(ch(j),'ZData');

  xdata(nan_mask) = nan;
  ydata(nan_mask) = nan;
  zdata(nan_mask) = nan;

  set(ch(j),'XData',xdata,'YData',ydata,'ZData',zdata);
end;
