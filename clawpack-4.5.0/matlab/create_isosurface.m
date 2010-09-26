function p = create_isosurface(xc,yc,zc,q,r,ivalue,icolor,ialpha,mappedgrid)

% Internal matlab routine used by Clawpack graphics.

[xcm,ycm,zcm] = meshgrid(xc,yc,zc);

p = patch(isosurface(xcm,ycm,zcm,q,ivalue));

reducepatch(p,r);

if (mappedgrid == 1)
  v = get(p,'Vertices');
  if (~isempty(v))
    [v(:,1),v(:,2),v(:,3)] = mapc2p(v(:,1),v(:,2),v(:,3));
    set(p,'Vertices',v);
    % isonormals are computed using triangular elements - results may not
    % appear smooth, but it isn't clear how to set isonormals in the mapped
    % grid case.
  end;
else
  % isonormals are computed using data - smoother results?
  isonormals(xcm,ycm,zcm,q,p);  % for lighting effects.
end;
if (strcmp(icolor,'q')==1)
  set(p,'CData',ivalue);
  set(p,'FaceColor','flat');
else
  [l,ic,m,msg] = colstyle(icolor);
  if (strcmp(ic,'') == 1)
    error(['*** create_isosurface : Bad color spec : ''',icolor,'''']);
  end;
  icolor = ic;
  set(p,'FaceColor',icolor);
end;
set(p,'EdgeColor','none');
set(p,'FaceAlpha',ialpha);
udata.q = q;
udata.xc = xc;
udata.yc = yc;
udata.zc = zc;
udata.value = ivalue;
udata.alpha = ialpha;
udata.color = icolor;
udata.mappedgrid = mappedgrid;
set(p,'UserData',udata);

set(p,'Tag','on');
