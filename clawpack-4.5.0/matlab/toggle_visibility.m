function toggle_visibility(p,obj_state,obj_handles)

% Internal matlab routine for Clawpack graphics.

% This sets obj_visibility, according to the following truth table:
%
% ---------------------------------------------------------
% obj_state  | patch_visibility_state  |  obj_visibility
% ---------------------------------------------------------
%     T      |         T               |       T
%     T      |         F               |       F
%     F      |         T               |       F
%     F      |         F               |       F
% ---------------------------------------------------------
%
% The end effect is that the graphics objects (patches, contourlines, patch
% borders,....) are only visibible if their patch is visible AND  the user
% defined state is 'on'.

if (isempty(obj_state) | isempty(obj_handles))
  return;
end;

% Get visibility state of current patch
p_status = strcmp(get(p,'Visible'),'on');

% Get value of truth table above, for array of values in obj_state.
set_visible = p_status & obj_state;

set(obj_handles(set_visible),'Visible','on');
set(obj_handles(~set_visible),'Visible','off');
