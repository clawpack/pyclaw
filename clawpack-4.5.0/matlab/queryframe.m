%  QUERYFRAME used by PLOTCLAW1, PLOTCLAW2 and PLOTCLAW3 to loop through data.
%
%      QUERYFRAME is normally called directly from PLOWCLAW<N>, and so is
%      not typically called by the user.  However, the user can supress the
%      query option by setting the 'NoQuery' parameter to 1.  This is used
%      for example, in making movies.
%
%      See also PLOTCLAW1, PLOTCLAW2, PLOTCLAW3.

if exist('NoQuery')
  if NoQuery == 1
    % set NoQuery=1 if you want the plots to be produced with no
    % query in between.  Particularly useful if you want lots of frames to
    % be printed out for an animation (put a command like makeframegif
    % in afterframe.m and set NoQuery=1)
    pause(1)
    Frame = Frame + 1;
    if Frame > MaxFrames
      break;   % break out of plotclawN after last frame
    end
    return
  end
end


inp = 'k';
while strcmp(inp,'k')

  inp = input(...
      ['Hit <return> for next plot, or type k, r, rr, j, i, q, or ?  '],'s');

  if strcmp(inp,'?')
    disp('  k  -- keyboard input.  Type any commands and then "return"')
    disp('  r  -- redraw current frame, without re-reading data')
    disp('  rr -- re-read current file,and redraw frame');
    disp('  j  -- jump to a particular frame')
    disp('  i  -- info about parameters and solution')
    disp('  q  -- quit')
  elseif strcmp(inp,'k')
    keyboard
  elseif strcmp(inp,'r')
    % redraw:  leave Frame counter alone
    if Frame==-1
      disp('Cannot redraw yet')
      inp = 'k';
    end
  elseif strcmp(inp,'rr')
    % redraw frame AND re-read data
    amrdata = [];
  elseif strcmp(inp,'j')
    Frame = input('Frame to jump to? ');
  elseif strcmp(inp,'i')
    if clawdim == 1
      infoplot1
      disp(' ')
      disp(' ')
      disp('hit <return> for information about this frame')
      pause
      infoframe1
    end
    if clawdim == 2
      infoplot2
      disp(' ')
      disp('hit <return> for information about this frame')
      pause
      infoframe2
    end
    if clawdim == 3
      infoplot3
      disp(' ');
      disp('hit <return> for information about this frame');
      pause
      infoframe3
    end;
    inp = 'k';
  elseif isempty(inp)
    % go to next frame
    Frame = Frame + 1;
  elseif (~strcmp(inp,'q'))
    % quit handled separately below.
    % Otherwise unrecognized input, go back and try again
    inp = 'k';
  end % if strcmp
end % while strcmp

if strcmp(inp,'q')
  % quit now
  break
end
