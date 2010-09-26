function h = getlinehandles(level)

% Internal matlab routine

if (nargin == 0)
  level = 1;
end;

if (length(level) > 1)
  error(['getlinehandles : You can only get handles for a single ',...
	'level of amr data at a time']);
end;

amrlines = get_lines;

h = amrlines{level};
