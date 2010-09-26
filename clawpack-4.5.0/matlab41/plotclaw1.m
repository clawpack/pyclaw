%
% plotclaw1.m
%
% generic plotting routine for claw1d and amrclaw1d output in matlab
% R. J. LeVeque, 1999
%
% Various parameters are set in setplot1.m
% The default version in claw/matlab/setplot1.m can be copied to your
% directory and modified to set things up differently, or type 'k'
% at the prompt to get keyboard control and change a value.
%
%---------------------------------------------------------------------

clawdim = 1;

disp(' ')
disp('plotclaw1  plots 1d results from clawpack')

% set plotting parameters:
whichfile = which('setplot1');
if strcmp(whichfile,'')
    disp('*** No setplot1 file found')
  else
    inp = input(['Execute setplot1 (default = no)? '],'s');
    if (strcmp(inp,'y') | strcmp(inp,'yes'))
       setplot1
       disp(' ')
       disp(['Executing m-script ' whichfile])
       end
  end
disp(' ')

% the file setprob.m can be used to set up any necessary physical parameters
% or desired values of plotting parameters for this particular problem.

whichfile = which('setprob');
if strcmp(whichfile,'')
    %disp('*** No setprob file found')
  else
    disp(['Executing m-script ' whichfile])
    disp(' ')
    setprob
  end


%=============================================
% MAIN LOOP ON FRAMES:
%=============================================

Frame = -1;  % initialize frame counter
if ~exist('MaxFrames')
   disp('MaxFrames parameter not set... you may need to execute setplot1')
   break
   end

while Frame <= MaxFrames

    % pause for input from user to determine if we go to next frame,
    % look at data, or skip around.  This may reset Frame counter.

    queryframe

    % produce the plot:

    plotframe1   % routine claw/matlab/plotframe1.m reads data and
                 % does the plotting for this frame

    end % main loop on frames
