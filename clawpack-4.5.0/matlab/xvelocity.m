function u = xvelocity(data)
%
% Compute the x-velocity from the data
% for problems where the second component is the x-momentum and the
% first component is the "density"  (e.g. Euler, shallow water equations)

u = data(:,2)./data(:,1);

