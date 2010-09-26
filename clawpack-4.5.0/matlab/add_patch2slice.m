function add_patch2slice(sdir,sval,snum,xc,yc,zc, xe,ye,ze, q,level,...
    contourlevels,mappedgrid,manifold,maskflag,grid_number);

% Internal matlab routine for Clawpack graphics.

% This routine adds a new patch to slice 'snum', in direction 'sdir', and sets
% various initial user defined settings.

slices = get_slices(sdir);
slice = slices{snum};

% Mask any patches that might be underneath the new patch we want to add.
if (maskflag == 1)
  mask_patches(slice,sdir,level,xe(1),xe(end),ye(1),ye(end),ze(1),ze(end));
end;

% Get the new patch
new_patch = create_patch(xc,yc,zc,xe,ye,ze,q,sdir,sval,contourlevels,...
    mappedgrid,manifold,grid_number);

% Finally add the patch to current slice, and update figure UserData.
slices{snum}{level}(end+1) = new_patch;
set_slices(sdir, slices);
