function set_isosurface_visibility(p,vstr);

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
