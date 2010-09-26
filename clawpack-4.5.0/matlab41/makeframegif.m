% Produce a gif file of this frame with numbering so that gifmerge
% can be used to make an animation:  frame00000.gif, frame00001.gif, etc.
%
% By putting the command makeframegif in the afterframe.m file,
% each frame will converted into a gif file when plotclaw1/2/3 is used.
% Then the unix command 
%    gifmerge frame*.gif > movie.gif
% should produce an animated gif file.

framest = num2str(Frame);
while(size(framest,2))<5,
  framest = ['0',framest];
  end;
fname = ['frame' framest];
printgif(fname);

