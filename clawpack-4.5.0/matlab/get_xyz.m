function varargout = get_xyz(varargin);

% Internal matlab routine for Clawpack graphics.

sdir = varargin{4};
idir = findstr(lower(sdir),'xyz');

if (idir == 1)
  iperm_rev = [1;2;3];
elseif (idir == 2)
  iperm_rev = [3;1;2];
else
  iperm_rev = [2;3;1];
end;

for i = 1:nargout,
  varargout{i} = varargin{iperm_rev(i)};
end;
