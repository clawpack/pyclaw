function set_cline_visibility(p)

% Internal matlab routine for Clawpack graphics.

if (~ishandle(p))
  return;
end;

p_status = get(p,'Visible');
udata = get(p,'UserData');
clines = udata.contourLines;

% Not every patch will have contourlines, even if ContourValues is
% specified, and patch is visible.
if (~isempty(clines))
  cline_state = strcmp(get(clines,'Tag'),'on');
  toggle_visibility(p,cline_state,clines);
end;
