function create_sliceintersections(mappedgrid)

% Internal matlab routine for Clawpack graphics.

% Now plot intersections of x-y planes, x-z planes and y-z planes.
sdirs = {'x', 'y', 'z'};
ve_names = {'xe','ye','ze'};
isect_names = {'yzIntersect','xzIntersect','xyIntersect'};

for idir = 1:3,
  for jdir = (idir+1):3,
    slices_idir = get_slices(sdirs{idir});
    slices_jdir = get_slices(sdirs{jdir});

    % Now consider all pair-wise combinations of these two sets of slices
    for n = 1:length(slices_idir);
      for m = 1:length(slices_jdir);

	% cell2mat doesn't seem to exist on some versions of Matlab
	% so I have copied my matlab version to ../matlab2 files.
	% pvec = cell2mat_claw(slices_idir{n});
	% qvec = cell2mat_claw(slices_jdir{m});

	pvec = [];
	qvec = [];
	% Need a big list of all patches we can intersect
	for ll = 1:length(slices_idir{n}),
	  pvec = [pvec slices_idir{n}{ll}];
	end;

	for ll = 1:length(slices_jdir{m}),
	  qvec = [qvec slices_jdir{m}{ll}];
	end;

	for k = 1:length(pvec),
	  for l = 1:length(qvec),
	    p = pvec(k);
	    q = qvec(l);

	    udata_p = get(p,'UserData');
	    udata_q = get(q,'UserData');

	    np = udata_p.grid_number;
	    nq = udata_q.grid_number;
	    if (np == nq)
	      % This sets two of the coordinates for the line we are going
	      % to draw
	      sval_p = udata_p.sval;
	      sval_q = udata_q.sval;

	      kdir = (1:3);
	      kdir([idir; jdir]) = [];

	      % Third direction is all that is left
	      ke = getfield(udata_p,ve_names{kdir});
	      vdata{idir} = sval_p + 0*ke;
	      vdata{jdir} = sval_q + 0*ke;
	      vdata{kdir} = ke;

	      if (mappedgrid == 1)
		[xdata, ydata, zdata] = mapc2p(vdata{:});
	      else
		[xdata, ydata, zdata] = deal(vdata{:});
	      end;

	      h = line('XData',xdata, 'YData',ydata,'ZData',zdata);

	      sp = getfield(udata_p,isect_names{kdir});
	      sp(end+1).line = h;
	      sp(end+1).sharedPatch = q;

	      sq = getfield(udata_q,isect_names{kdir});
	      sq(end+1).line = h;
	      sq(end+1).sharedPatch = p;

	      udata_p = setfield(udata_p,isect_names{kdir},sp);
	      udata_q = setfield(udata_q,isect_names{kdir},sq);
	      set(p,'UserData',udata_p);
	      set(q,'UserData',udata_q);
	    end;  % np == nq
	  end; % (l) loop over all patches
	end;   % (k) loop over all patches
      end; % (m) loop over all slices
    end;   % (n) loop over all slices
  end;  % jdir
end; % idir
