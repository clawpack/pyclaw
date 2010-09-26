function ...
    intersectQWithSlices(xg,yg,zg,xe,ye,ze,q,slice_dir,slice_list,...
    level,plotgrid,plotgridedges,plottype,contour_levels);

global patch_handles patch_count;
global contour_handles;

xlow  = xe(1);
xhigh = xe(end);
ylow  = ye(1);
yhigh = ye(end);
zlow  = ze(1);
zhigh = ze(end);

% Loop over constant slices.
for slice_i = 1:length(slice_list),

  slice_const = slice_list(slice_i);  % constant value (either x,y or z) of slice.

  % Check to see that 3D patch intersects with slice_i
  if (strcmp(slice_dir,'x'))
    isect = xlow <= slice_const & slice_const <= xhigh;
  elseif (strcmp(slice_dir,'y'))
    isect = ylow <= slice_const & slice_const <= yhigh;
  elseif (strcmp(slice_dir,'z'))
    isect = zlow <= slice_const & slice_const <= zhigh;
  end;

  if (isect == 1)
    patch_count(slice_i,level) = patch_count(slice_i,level) + 1;
    if (level > 1)
      patches_underneath = patch_handles(slice_i,1:(level-1),:);
      lines_underneath = contour_handles(slice_i,1:(level-1),:);
      maskPatches(xlow,xhigh,ylow,yhigh,zlow,zhigh,...
	  patches_underneath,lines_underneath,~isempty(contour_levels),slice_dir);
    end;
    [new_patch,new_lines] = plotSlicePatch(xg,yg,zg,xe,ye,ze,q,slice_dir,...
	slice_const,plotgrid,contour_levels);

    if (plottype == 2)
      set(new_patch,'FaceColor','w');
    end;

    if (plotgridedges)
      h = projectBorder(xlow,xhigh,ylow,yhigh,zlow,zhigh,slice_dir,slice_const);
      userdata = get(new_patch,'UserData');
      userdata = setfield(userdata,'Border',h);
      set(new_patch,'UserData',userdata);
    end;

    patch_num = patch_count(slice_i,level);
    patch_handles(slice_i,level,patch_num) = new_patch;
    contour_handles{slice_i,level,patch_num} = new_lines;

  end; % end intersection of cube with slice
end; % end for loop over slices.
