function [linestyle,linecolors,markerstyle] = get_plotstyle(pstyle,len)

% Internal Matlab routine for Clawpack graphics

lstyle = {};
mstyle = {};
lcolors = {};
for i = 1:length(pstyle),
  [l,c,m,msg] = colstyle(pstyle{i});
  if (~isempty(msg))
    error(msg);
  end;
  if (isempty(l))
    lstyle{i} = 'none';
  else
    lstyle{i} = l;
  end;
  if (isempty(m))
    mstyle{i} = 'none';
  else
    mstyle{i} = m;
  end;
  if (isempty(c))
    lcolors{end+1} = 'b';
  else
    lcolors{end+1} = c;
  end;
end;

linestyle = set_length(lstyle,len);
markerstyle = set_length(mstyle,len);
linecolors = set_length(lcolors,len);
