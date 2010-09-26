function h = create_clines(c,sval,sdir,mappedgrid,manifold)

% Internal matlab routine for Clawpack graphics.


% This draws contour lines stored in contour matrix c.  THis may be
% called by either 2d or 3d routines.

global creatingclines

if (isempty(c))
  % No contour lines on this patch
  h = [];
  return;
end;

st_idx = 1;
line_cnt = 0;
while (1)
  cval = c(1,st_idx);
  next_length = c(2,st_idx);
  line_cnt = line_cnt + 1;
  x_like = zeros(1,next_length) + sval;
  y_like = c(1,st_idx+1:st_idx+next_length);
  z_like = c(2,st_idx+1:st_idx+next_length);

  [xdata,ydata,zdata] = get_xyz(x_like,y_like,z_like,sdir);

  udata.cartCoords = [xdata', ydata', zdata'];
  if (mappedgrid == 1 | manifold == 1)
    if (mappedgrid == 1)
      if (nargin('mapc2p') == 2)
	[xdata,ydata] = mapc2p(xdata,ydata);
      else
	[xdata, ydata, zdata] = mapc2p(xdata,ydata,zdata);
      end;
    end;
    if (manifold == 1)
      creatingclines = 1;  % flag for mapc2m to indicate contour lines
                           % are being generated: may want to shift off
                           % manifold slightly to avoid hidden line removal
                           % causing apparent gaps in contours.
      [xdata, ydata, zdata] = mapc2m(xdata,ydata);
      creatingclines = 0;
    end;
  end;

  h(line_cnt) = line('XData',xdata,'YData',ydata,'ZData',zdata,'Color','k');
  set(h(line_cnt),'Tag','on');
  udata.cval = cval;
  set(h(line_cnt),'UserData',udata);
  st_idx = st_idx + next_length + 1;
  if (st_idx > length(c))
    return;
  end;
end;
