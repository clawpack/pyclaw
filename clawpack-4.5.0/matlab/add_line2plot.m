function add_line2plot(s,q,level,mstyle,lcolor,lstyle)

% Internal Matlab routine for Clawpack graphics.

amrlines = get_lines;

h = line(s,q,'Marker',mstyle,'Color',lcolor,'LineStyle',lstyle);

udata.s = s;
udata.q = q;
set(h,'UserData',udata);

amrlines{level}(end+1) = h;
set_lines(amrlines);
