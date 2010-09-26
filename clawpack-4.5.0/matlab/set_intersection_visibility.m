function set_intersection_visibility(p);

% Internal matlab routine for Clawpack graphics.

if (~ishandle(p))
  return;
end;

udata = get(p,'UserData');
names = {'xyIntersect','xzIntersect','yzIntersect'};
for i = 1:3,
  f = getfield(udata,names{i});
  if (~isempty(f))
    line_state = strcmp(get([f.sharedPatch],'Visible'),'on');
    lines = [f.line];
    toggle_visibility(p,line_state,lines);
  end;
end;
