function set_patch_visibility(p,vstr);

% Internal matlab routine for Clawpack graphics.

if (~ishandle(p))
  return;
end;

patch_state = strcmp(get(p,'Tag'),'on');
if (patch_state == 1)
  set(p,'Visible',vstr);
else
  set(p,'Visible','off');
end;

% We are toggling the visibility status from on-->off or off-->on
% set(p,'Visible',vstr);

% This will use visibility of current patch to determine whether to turn
% graphics objects on or off.

set_cline_visibility(p);   % contourlines
set_pborder_visibility(p); % patch borders
set_mesh_visibility(p);    % coarsened mesh
set_intersection_visibility(p);   % 3d only;  for intersection lines between slices.
