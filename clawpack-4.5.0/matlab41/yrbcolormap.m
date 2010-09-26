%
% yellow-red-blue color map 
%
c01 = 0:.1:1;
c10 = 1:-.1:0;
c0 = 0*c01;
c1 = 1 + c0;
yrb = [c1' c10' c0'; c10' c0' c01'];
colormap(yrb)

