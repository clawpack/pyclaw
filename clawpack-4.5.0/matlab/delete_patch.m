function delete_patch(p)

% Internal matlab routine for Clawpack graphics.

if (~ishandle(p))
  return;
end;

udata = get(p,'UserData');

if (ishandle(udata.contourLines))
  delete(udata.contourLines);
end;

if (ishandle(udata.border))
  delete(udata.border);
end;

if (ishandle(udata.mesh.xlines))
  delete(udata.mesh.xlines);
end;

if (ishandle(udata.mesh.ylines))
  delete(udata.mesh.ylines);
end;

% Now get the rest lines intersecting slices...
names = {'xyIntersect','xzIntersect','yzIntersect'};
for i = 1:3,
  xyzlines =  getfield(udata,names{i});  % dynamic structure field
  for n = 1:length(xyzlines),
    if (ishandle(xyzlines(n).line))
      % We can only delete this handle if no other patch has already
      % deleted it.
      delete(xyzlines(n).line); % Get all lines associated with this patch
    end;
  end;
end;
delete(p);  % Delete current patch
