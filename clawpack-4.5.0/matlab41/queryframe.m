% used by plotclaw1.m and plotclaw2.m to query user about what to do next

    inp = 'k';
    while strcmp(inp,'k')

      inp = input(...
	['Hit <return> for next plot, or type k, r, j, i, q, or ?  '],'s');

      if strcmp(inp,'?')
         disp('  k -- keyboard input.  Type any commands and then "return"')
         disp('  r -- redraw current frame')
         disp('  j -- jump to a particular frame')
         disp('  i -- info about parameters and solution')
         disp('  q -- quit')
      elseif strcmp(inp,'k')
         keyboard
      elseif strcmp(inp,'r')
	 % redraw:  leave Frame counter alone
	 if Frame==-1
	    disp('Cannot redraw yet')
	    inp = 'k';
	    end
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

