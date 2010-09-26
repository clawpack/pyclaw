function set_pborder_visibility(p)

% Internal matlab routine for Clawpack graphics.

if (~ishandle(p))
  return;
end;

p_status = get(p,'Visible');
udata = get(p,'UserData');

border_state = strcmp(get(udata.border,'Tag'),'on');
toggle_visibility(p,border_state,udata.border);
