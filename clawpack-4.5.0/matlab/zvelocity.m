function w = zvelocity(data)
%
% Compute the z-velocity from the data
% for problems where the fourth component is the z-momentum and the
% first component is the "density"  (e.g. Euler)


w = data(:,4)./data(:,1);
