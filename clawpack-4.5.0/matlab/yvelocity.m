function v = yvelocity(data)
%
% Compute the y-velocity from the data
% for problems where the third component is the y-momentum and the
% first component is the "density"  (e.g. Euler, shallow water equations)

v = data(:,3)./data(:,1);

