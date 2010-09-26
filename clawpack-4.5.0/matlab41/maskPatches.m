function maskPatches(xlow,xhigh,ylow,yhigh,zlow,zhigh,punder,cunder,...
    clines,sdir,sval)

% This masks out any patches that are directly underneath the current
% set of patches to be plotted.  Usually, these will only be patches on the
% next level down (level = 'current level' - 1), but if data from that level is
% not being plotted for some reason (perhaps level 1 is not plotted, if it
% is too time consuming) then we must go to level = 'current level' - 2, etc.

% clines == 1 if we have contour lines; otherwise it is zero.

% Dimensions of punder: For this particular slice, we have 'nlevels' levels,
% and at each level, we have 'npatches' patches.  'nlevels' should be equal
% to 'current level' - 1.


[nlevels,npatches] = size(punder);

mask_level = nlevels;
% Find first non-zero level directly beneath current level.
while (punder(mask_level) == 0) % not a patch handle
  mask_level = mask_level-1;
  if (mask_level == 0)
    % There are no patches under this one to mask
    return;
  end;
end;

patch_idx = find(punder(mask_level,:) ~= 0);

for i = 1:length(patch_idx)
  patch_i = patch_idx(i);
  p = punder(mask_level,patch_i);

  % List of vertices in patch p
  v = get(p,'Vertices');

  xv = v(:,1); % x values
  yv = v(:,2); % y values
  zv = v(:,3); % z values

  % Mask for bad vertices

  % Make sure that masked region really is in mask.
  s = 1e-5;
  if (strcmp(sdir,'x') == 1)
    m = (yv > ylow+s & yv < yhigh-s & zv > zlow+s & zv < zhigh-s);
  elseif (strcmp(sdir,'y') == 1)
    m = (xv > xlow+s & xv < xhigh-s & zv > zlow+s & zv < zhigh-s);
  elseif (strcmp(sdir,'z')== 1)
    m = (xv > xlow+s & xv < xhigh-s & yv > ylow+s & yv < yhigh-s);
  end;

  % indices of bad vertices
  masked_vertices = find(m);
  v(masked_vertices,:) = nan;

  % Reset vertices in patch to new values, with NaN's.  These won't get
  % plotted on next redraw of patch.
  set(p,'Vertices',v);

  % Now mask any contour lines
  if (clines == 1)
  % if (strcmp(get(cunder{mask_level,patch_i},'Type'),'line')==1)
    ch = cunder{mask_level,patch_i};
    for j = 1:length(ch),  % loop over contour lines on this patch
      xdata = get(ch(j),'XData');
      ydata = get(ch(j),'YData');
      zdata = get(ch(j),'ZData');
      if (strcmp(sdir,'x') == 1)
	m = (ydata > ylow+s & ydata < yhigh-s & zdata > zlow+s & zdata < zhigh-s);
      elseif (strcmp(sdir,'y') == 1)
	m = (xdata > xlow+s & xdata < xhigh-s & zdata > zlow+s & zdata < zhigh-s);
      elseif (strcmp(sdir,'z')== 1)
	m = (xdata > xlow+s & xdata < xhigh-s & ydata > ylow+s & ydata < yhigh-s);
      end;
      xdata(find(m)) = nan;
      ydata(find(m)) = nan;
      zdata(find(m)) = nan;
      set(ch(j),'XData',xdata,'YData',ydata,'ZData',zdata);
    end;
  end;

end;  % End loop on patches.
