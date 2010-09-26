function showmesh(ratio,level,sdir,snum)

% SHOWMESH shows a coarsened mesh on specified amr level.
%
%      SHOWMESH(RATIO) shows a grid coarsened by a factor of RATIO
%      relative to given level, on all amr levels.
%
%      SHOWMESH(RATIO,LEVEL) shows a coarsened mesh on levels specified in
%      the vector LEVEL.
%
%      SHOWMESH, by itself, shows the computational grid on all levels, and is
%      has the same effect as SHOWGRIDLINES().
%
%      For 2d manifold plots, or 3d slice plots, the mesh on a give level
%      is not masked by amr patches on higher levels, so for example,
%      SHOWMESH(4,1) has the effect of showing a mesh coarsened by factor
%      of 4, relative to level 1, over the entire manifold or slice. On 2d
%      flat plots, the mesh lines are masked.
%
%      See also HIDEMESH, SHOWGRIDLINES, HIDEGRIDLINES.

if (nargin < 3)
  sdirs = {'x','y','z'};
  if (nargin < 2)
    level = 1;
    if (nargin < 1)
      ratio = 1;
    end;
  end;
  hidemesh;  % hide all levels
else
  sdirs = {sdir};
  hidemesh(level);
end;

for idir = 1:length(sdirs),
  slices = get_slices(sdirs{idir});
  if (nargin < 4)
    snum = 1:length(slices);
  end;
  for ns = 1:length(snum),
    n = snum(ns);
    if (n < 1 | n > length(slices))
      continue;
    end;
    slice = slices{n};
    for l = 1:length(level),
      pvec = slice{level(l)};
      for k = 1:length(pvec),
	p = pvec(k);
	set_blocknumber(k);  % In case we are doing a multiblock plot
	udata = get(p,'UserData');
	if (isempty(udata.mesh.xlines) | isempty(udata.mesh.ylines))
	  line_dir = 1:2;  % do both x and y lines.
	  mesh = ...
	      create_mesh(sdirs{idir},udata.sval,line_dir, ...
	      udata.xe,udata.ye,udata.ze,...
	      udata.mappedgrid,udata.manifold);
	  udata.mesh = mesh;
	  set(p,'UserData',udata);
	else
	  mesh = udata.mesh;
	end;
	if (nargin >= 1)
	  udata.mesh.ratio = ratio;
	  set(p,'UserData',udata);
	end;
	xlines = mesh.xlines;
	ylines = mesh.ylines;
	lx = length(xlines);
	ly = length(ylines);

	xmask = kron(ones(lx,1),[1;zeros(ratio-1,1)])==1;
	xmask(lx) = 1==1;  % 'true' not in earlier versions of Matlab
	set(xlines( xmask(1:lx)),'Tag','on');
	set(xlines(~xmask(1:lx)),'Tag','off');

	ymask = kron(ones(ly,1),[1;zeros(ratio-1,1)])==1;
	ymask(ly) = 1==1;  % 'true' not earlier versions of Matlab
	set(ylines( ymask(1:ly)),'Tag','on');
	set(ylines(~ymask(1:ly)),'Tag','off');

	set_mesh_visibility(pvec(k));
      end;
    end;
  end;
end;
