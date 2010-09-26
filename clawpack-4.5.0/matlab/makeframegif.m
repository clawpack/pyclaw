% MAKERAMEGIF produces a .gif file of the current frame.
%
%     Produce a gif file of this frame with numbering so that gifmerge
%     can be used to make an animation:  frame00000.gif, frame00001.gif, etc.
%
%     By putting the command makeframegif in the afterframe.m file,
%
%     Then catenate all these files together into an animation.
%     There are various ways to do this in unix or linux depending on what
%     software is on your system, e.g.
%         convert -delay 20 frame*.gif movie.gif
%         gifmerge frame*.gif > movie.gif
%         gifsicle frame*.gif > movie.gif
%     gifsicle for linux can be freely obtained from
%             http://www.lcdf.org/gifsicle/
%
%     See also PRINTGIF, PRINTJPG.
%

framest = num2str(Frame);
while(size(framest,2))<5,
  framest = ['0',framest];
  end;
fname = ['frame' framest];
printgif(fname);
