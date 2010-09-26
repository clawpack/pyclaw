function set_mesh_visibility(p);

% Internal matlab routine for Clawpack graphics.

if (~ishandle(p))
  return;
end;

udata = get(p,'UserData');
xlines = udata.mesh.xlines;
xline_state = strcmp(get(xlines,'Tag'),'on');
toggle_visibility(p,xline_state,xlines);

ylines = udata.mesh.ylines;
yline_state = strcmp(get(ylines,'Tag'),'on');
toggle_visibility(p,yline_state,ylines);
