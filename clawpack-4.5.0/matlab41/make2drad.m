function q2 = make2drad(x0,y0,r1,q1,x2,y2)
%
% function q2 = make2drad(x0,y0,r1,q1,x2,y2)
%
% make2drad:  Convert 1d data into a radially-symmetric
%             2d grid function by interpolation
%
%  on input:  r1 is a vector of radii with r1(1)=0
%             q1 is a vector of values q(r1)
%             (x0,y0) is point about which 2d solution should be centered
%             x2 is array of x values on 2d grid
%             y2 is array of y values on 2d grid
%
%  on output: q2(i,j) = q evaluated at
%                 r(i,j) = sqrt((x2(i,j)-x0)^2 + y2(i,j)-y0)^2)

if size(x2) ~= size(y2) 
   error('Error ****  need size(x2)=size(y2)')
   end

m1 = length(r1);

r2 = sqrt((x2-x0).^2 + (y2-y0).^2);
for i=1:size(x2,1)
   for j=1:size(x2,2)
      q2(i,j) = interp1(r1,q1,r2(i,j));
      q2(isnan(q2)) = q1(m1);
      end % loop on j
   end % loop on i

