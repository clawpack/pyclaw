%
% red-yellow--blue color map 
%
c01 = 0:.1:1;
c10 = 1:-.1:0;
c0 = 0*c01;
c1 = 1 + c0;
mymap = [c1' c0' c0'; c1' c01' c0'; c10' c10' c01';   c0' c0' c1'];
colormap(mymap)
