%
% If you have a file called aftergrid.m on your matlab path, the commands
% in this file will be executed after each grid of data is plotted when
% plotting AMR results.  This is useful if you want to pause after plotting
% each grid to manipulate the data (put "keyboard" command in aftergrid.m)
% or if you want to label each grid with it's gridno, using for example
%   text(mean(xcenter),mean(ycenter),num2str(gridno))
