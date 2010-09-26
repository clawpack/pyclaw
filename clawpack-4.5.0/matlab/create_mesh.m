function mesh = create_mesh(sdir,sval,line_dir,xe,ye,ze,mappedgrid,manifold)

% Internal matlab routine for Clawpack graphics.

% This routine is SLOW because I draw all the lines, regardless of whether
% they are going to be shown by the showmesh command.  A better way to do
% this would be to return a structure which stores the (xdata,ydata,zdata)
% that would be used if the lines were drawn;  Then let the showmesh command
% selectively draw lines which should be turned on and off.  But this might
% be hard...

[xe_like, ye_like, ze_like] = get_xyzlike(xe,ye,ze,sdir);

mesh.xlines = [];
mesh.ylines = [];
ve = {ye_like, ze_like};
names = {'xlines','ylines'};

for m = 1:length(line_dir)
  n = line_dir(m);
  np1 = mod(n,2) + 1;
  vlines = [];
  for k = 1:length(ve{n}),
    yze_like = {0*ve{np1} + ve{n}(k) ve{np1}};
    xe_like = 0*ve{np1} + sval;
    ye_like = yze_like{n};
    ze_like = yze_like{np1};

    [xdata, ydata, zdata] = get_xyz(xe_like,ye_like,ze_like,sdir);

    udata.cartCoords = [xdata' ydata', zdata'];
    if (mappedgrid == 1 | manifold == 1)
      if (mappedgrid == 1)
	if (nargin('mapc2p') == 2)
	  [xdata,ydata] = mapc2p(xdata,ydata);
	else
	  [xdata, ydata, zdata] = mapc2p(xdata,ydata,zdata);
	end;
      end;
      if (manifold == 1)
	[xdata, ydata, zdata] = mapc2m(xdata,ydata);
      end;
    end;
    hdl = line('XData',xdata,'YData',ydata,'ZData',zdata);
    set(hdl,'Visible','off');
    set(hdl,'Tag','off');
    set(hdl,'UserData',udata);
    vlines(k) = hdl;
  end;
  mesh = setfield(mesh,names{n},vlines);
end;
mesh.ratio = 1;
