function h = plotCube(x1,x2,y1,y2,z1,z2)

% This plots a box in 1, 2, or 3d.

h(1) = line('XData',[x1 x1 x2 x2 x1],'YData',[y1 y2 y2 y1 y1],...
    'ZData',ones(1,5)*z1,'Color','k');

h(2) = line('XData',[x1 x1 x2 x2 x1],'YData',[y1 y2 y2 y1 y1],...
    'ZData',ones(1,5)*z2,'Color','k');

h(3) = line('XData',ones(1,5)*x1,'YData',[y1 y2 y2 y1 y1],...
    'ZData',[z1 z1 z2 z2 z1],'Color','k');

h(4) = line('XData',ones(1,5)*x2,'YData',[y1 y2 y2 y1 y1],...
    'ZData',[z1 z1 z2 z2 z1],'Color','k');
