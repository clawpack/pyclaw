function [xm,ym,zm] = mapc2m(xc,yc)
% maps to points on a sphere of radius r

ijlower = find(xc<-1);  %indices of points on lower hemisphere
xc(ijlower) = -2-xc(ijlower); %flip across the line x=-1

% compute xp and yp by mapping [-1,1]x[-1,1] to unit circle
xc1 = abs(xc);
yc1 = abs(yc);
d = max(xc1,yc1);
d = max(d,1.e-10);
r1 = 1.; 
 
D = r1*d.*(2*ones(size(d))-d)./sqrt(2);

R = r1*ones(size(d));
center = D-sqrt(R.^2. - D.^2.);
xp = D.*xc1./d;
yp = D.*yc1./d;
ij = find(yc1==d);
yp(ij)=center(ij)+sqrt(R(ij).^2. - xp(ij).^2.);
ij = find(xc1==d);
xp(ij)=center(ij)+sqrt(R(ij).^2. - yp(ij).^2.);
xp = sign(xc).*xp;
yp = sign(yc).*yp; 


zp = sqrt(1.-(xp.^2.+yp.^2.));
zp(ijlower) = -zp(ijlower);  %negate z in lower hemisphere 

xm = xp;
ym = yp;
zm = zp;
