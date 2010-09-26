% MAKERAMEJPG produces a .jpg file of the current frame.
%
%     Produce a jpg file of this frame with numbering so that they
%     can be used to make an animation:  frame00000.jpg, frame00001.jpg, etc.
%
%     By putting the command makeframejpg in the afterframe.m file,
%
%     Then catenate all these files together into an animation.
%     There are various ways to do this in unix or linux depending on what
%     software is on your system, e.g.
%         convert -delay 10 -adjoin frame*.jpg movie.mpeg
%         convert -resize 600x600 -delay 10 frame*.jpg movie.gif
%

framest = num2str(Frame);
while(size(framest,2))<5,
  framest = ['0',framest];
  end;
fname = ['frame' framest];
printjpg(fname);
