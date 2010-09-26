function mask_patches(slice,sdir,level, xlow,xhigh,ylow,yhigh,zlow,zhigh)

% Internal matlab routine for Clawpack graphics.

% This masks out any patches that are directly underneath the current
% set of patches to be plotted.

if (level == 1)
  % Nothing to mask in this case
  return;
end;

% Find first non-zero level directly beneath current level.
mask_level = level - 1;
while (isempty(slice{mask_level})) % no patches have been plotted at this level
  mask_level = mask_level-1;
  if (mask_level == 0)
    % There are no patches under this one to mask
    return;
  end;
end;

% Vector of patch handles to mask at level 'mask_level'.
patches_to_mask = slice{mask_level};

for k = 1:length(patches_to_mask)
  p = patches_to_mask(k);

  mask_patch(p,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh);
  mask_clines(p,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh);
  % mask_mesh(p,sdir,xlow,xhigh,ylow,yhigh,zlow,zhigh);

end;  % End loop on patches that need to be masked.
