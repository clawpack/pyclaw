function aout = set_length(a,m)

% Internal Matlab routine for Clawpack graphics

% This routine expands the length of a to a vector of length m.
% If length(a) >= m, aout(1:m) = a(1:m).
% If length(a) < m, then aout cycles through entries of a to fill up
% m entries.
%

if (size(a,1) ~= 1 & size(a,2) ~= 1)
  error('set_length : At least one dimension must be 1');
end;

n = length(a);
for i = 1:m,
  aout(i) = a(mod(i-1,n) + 1);
end;
