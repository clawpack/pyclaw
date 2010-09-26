function printjpg(fname)

% PRINTJPG prints a figure as a .jpg
%
%     PRINTJPG(FNAME) makes a small 3in x 3in figure suitable for
%     putting on a webpage, and that looks better than what is obtained
%     by shrinking down the standard output from print -djpg
%
% See also MAKEFRAMEGIF, PRINTJPG, and the unix command CONVERT.

set(gcf,'paperunits','inches','paperposition',[0 0 3 3])
eval(['print -djpeg90 ' fname]);
