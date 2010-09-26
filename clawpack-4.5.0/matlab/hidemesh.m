function hidemesh(level)

% HIDEMESH hides a coarsened mesh on specified amr level.
%
%      HIDEMESH(LEVEL) hides the coarsened mesh on levels specified in
%      the vector LEVEL.
%
%      HIDEMESH, by itself, hides the mesh shown with SHOWMESH on all levels
%
%      See also SHOWMESH, SHOWGRIDLINES, HIDEGRIDLINES.

sdir = {'x', 'y', 'z'};

for idir = 1:3,
  slices = get_slices(sdir{idir});
  for n = 1:length(slices),
    slice = slices{n};
    if (nargin < 1)
      level = 1:length(slice);
    end;
    for l = 1:length(level),
      pvec = slice{level(l)};
      for k = 1:length(pvec),
	p = pvec(k);
	udata = get(p,'UserData');
	if (~isempty(udata.mesh))
	  set(udata.mesh.xlines,'Tag','off');
	  set(udata.mesh.ylines,'Tag','off');
	  set_mesh_visibility(pvec(k));
	end;
      end;
    end;
  end;
end;
