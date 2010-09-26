function add_isosurface2plot(xc,yc,zc,q,level,...
    isosurfvalues,isosurfcolors,isosurfalphas,mappedgrid)

% Internal Matlab routine used by Clawpack graphics

isurf = get_isosurfaces;

q = permute(q,[2 1 3]);

for n = 1:length(isosurfvalues),
  ivalue = isosurfvalues(n);
  icolor = isosurfcolors(n);
  ialpha = isosurfalphas(n);
  p = create_isosurface(xc,yc,zc,q,1,ivalue,icolor,ialpha,mappedgrid);
  isurf{n}{level}(end+1) = p;
end;

set_isosurfaces(isurf);
