function surfloop()

% SURFLOOP loops over user-specified isosurfaces.
%
%     SURFLOOP loops over isosurfaces specified in vector IsosurfValues.
%
%     At start of loop, all isosurfaces are hidden.  The user then steps through
%     all isosurfaces one at a time. .
%
%     See also SHOWSURFS, HIDESURFS.


isurfaces = get_isosurfaces;

if (length(isurfaces) == 0)
  fprintf('isoloop : %sIsurfValues == []\n',sdir);
  return;  % Nothing to loop over
end;

% First hide all slices in direction dir.
hidesurfs;

notdone = 1;
next_surf = 0;
last_surf = 0;
while (notdone)
  s = input('Hit <return> for next surface, or type k, r, j, i, q, or ? ','s');

  if (isempty(s))
    next_surf = mod(next_surf+1,length(isurfaces)+1);
  elseif (strcmp(s,'j'))
    next_surf = input('Input surf number : ');
  elseif (strcmp(s,'k'))
    keyboard;
  elseif (strcmp(s,'q'))
    return;
  end;
  hidesurfs(last_surf);
  fprintf('\n');
  fprintf('Showing surf %d\n',next_surf);
  showsurfs(next_surf);

  last_surf = next_surf;

end;
