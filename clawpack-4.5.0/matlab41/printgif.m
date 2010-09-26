function printgif(fname)
% Print figure as a gif file fname.gif
% First prints as a ppm file and then converts using unix commands
% ppmquant and ppmtogif
%
% After creating many such gif files, e.g. frame00001.gif, frame00002.gif,  etc.
% these can be merged together to make an animation using the unix command
%   gifmerge frame*.gif > movie.gif
%
% See also makeframegif.m  for numbering by frames 

% temporary file name for first ppm output:
fname1 = [fname 'TMP'];  
str = ['print ' fname1 ' -dppmraw'];
eval(str);

% reduce to 256 colors or ppmtogif might croak:
str = ['ppmquant 256 ' fname1 '.ppm > ' fname '.ppm' ];
unix(str);

% convert to gif file:
str = ['ppmtogif ' fname '.ppm  > ' fname '.gif'];
unix(str);

% clean up by removing ppm files:
str = ['rm ' fname '.ppm '];
unix(str);
str = ['rm ' fname1 '.ppm '];
unix(str);

