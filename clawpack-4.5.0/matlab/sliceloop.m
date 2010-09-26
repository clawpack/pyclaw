function sliceloop(sdir)

% SLICELOOP loops over user-specified 3d slices in given direction.
%
%     SLICELOOP(SDIR) loops over slices in direction specified by SDIR.
%     SDIR may be 'x', 'y', or 'z'.
%
%     SLICELOOP, by itself, prompts the user for a direction (either x,y or z)
%     and loops over slices in that direction.
%
%     At start of loop, all slices in specified direction are hidden.
%     Then slices are viewed, one at a time, in specified direction.
%     The visibility of slices in any of the other two directions
%     will not be affected.
%
%     SLICELOOP only loops over those slices specified by variables
%     xSliceCoords, ySliceCoords or zSliceCoords, in setplot3.m.
%
%     See also SHOWSLICES, HIDESLICES, SETPLOT3.


if (nargin == 0)
  sdir = input(['Enter x, y, or z to specify direction over which to',...
      ' loop : '],'s');
end;

setviews;

slices = get_slices(sdir);

if (length(slices) == 0)
  fprintf('sliceLoop : %sSliceCoords == []\n',sdir);
  return;  % Nothing to loop over
end;

% First hide all slices in direction dir.
hideslices(sdir);

notdone = 1;
next_slice = 0;
last_slice = 0;
while (notdone)
  s = input('Hit <return> for next slice, or type k, r, j, i, q, or ? ','s');

  if (isempty(s))
    next_slice = mod(next_slice+1,length(slices)+1);
  elseif (strcmp(s,'j'))
    next_slice = input('Input slice number : ');
  elseif (strcmp(s,'k'))
    keyboard;
  elseif (strcmp(s,'q'))
    return;
  end;
  hideslices(sdir,last_slice);
  fprintf('\n');
  fprintf('Showing slice %d\n',next_slice);
  showslices(sdir,next_slice);

  last_slice = next_slice;
  if exist('afterslice')==2
    % make an m-file with this name for any other commands you
    % want executed at the end of drawing each slice, for example
    % to print a gif file for use in making an animation of a sweep-through
    afterslice;
    end

end;
