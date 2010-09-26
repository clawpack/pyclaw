function varargout = get_xyzlike(varargin);

% Internal matlab routine for Clawpack graphics.

sdir = varargin{4};
idir = findstr(lower(sdir),'xyz');
% iperm = circshift_claw((1:3)',-(idir - 1));
if (idir == 1)
  iperm = [1;2;3];
elseif (idir == 2)
  iperm = [2;3;1];
else
  iperm = [3;1;2];
end;

for i = 1:nargout,
  varargout{i} = varargin{iperm(i)};
end;
