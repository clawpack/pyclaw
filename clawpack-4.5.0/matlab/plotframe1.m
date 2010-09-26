%
% PLOTFRAME1 plots data from a single Clawpack output file.
%
%    PLOTFRAME1 is called from PLOTCLAW1, the driver script for Clawpack
%    graphics.  PLOTFRAME1 uses parameters, defined in the user's
%    workspace, usually in the SETPLOT1 script, to determine what kind of
%    plot to create, and various other plotting options.
%
%    See SETPLOT for a complete list of the parameters that the user can
%    set.
%
%    The basic plot type available in 1d is a 1d line plot, with user
%    specified symbols.
%
%    By specifying MAPPEDGRID == 1, a user defined function mapc2p will be
%    called
%
%    See SETPLOT, PLOTCLAW1, MAPPEDGRID.



if length(amrdata) == 0
  disp(' ');
  disp(['Frame ',num2str(Frame),' does not exist ***']);
  return;
end

disp(' ')
disp(['Frame ',num2str(Frame),' at time t = ',num2str(t)]);


if exist('beforeframe')==2
  beforeframe  % make an m-file with this name for any other commands you
  % want executed before drawing each frame, for example
  % if you want to use axes to specify exactly where the
  % plot will be in the window, aspect ratio, etc.
end

set_value('maxlevels','MaxLevels',6);

if (~exist('PlotStyle') & exist('plotstyle'))
  % Parse line spec style
  disp([' *** plotframe1 : ''plotstyle'' should be replaced by ',...
        '''PlotStyle''.  Set SETPLOTSTYLE.']);
  if (~iscell(plotstyle))
    % This provides backwards compatibility for plots with a single type.
    PlotStyle = {plotstyle};
  else
    PlotStyle = plotstyle;
  end;
end;

set_value('pstyle','PlotStyle',{'b-'});
if (~iscell(PlotStyle))
  error([' *** plotframe1 : PlotStyle must be a cell matrix. ',...
	'Use ''setplotstyle''.']);
end;

[linestyle,linecolors,markerstyle] = get_plotstyle(pstyle,maxlevels);

set_value('mappedgrid','MappedGrid',0);
if (mappedgrid == 1)
  if (~exist('mapc2p'))
    error(' *** plotframe1 : MappedGrid == 1 but no mapc2p function was found');
  end;
end;

clear_amrplot;
create_amrplot(maxlevels);

% Call to 'newplot' is needed in case the user has added any extra lines
% to the plot (in afterframe, for example), using 'hold on/off'.  'newplot'
% respects the hold on/off commands, whereas clear_amrplot and
% create_amrplot do not.
newplot;


%=============================================
% MAIN LOOP ON GRIDS FOR THIS FRAME:
%=============================================

qmin = [];
qmax = [];
for ng = 1:length(amrdata),

  gridno = amrdata(ng).gridno;
  level = amrdata(ng).level;
  mx = amrdata(ng).mx;
  xlow = amrdata(ng).xlow;
  dx = amrdata(ng).dx;

  data = amrdata(ng).data';

  if UserVariable==1
    % User has supplied a function to convert original q variables to
    % the variable which is to be plotted, e.g. Mach number, entropy.
    q = feval(UserVariableFile,data);
  else
    q = data(:,mq);
  end;

  amrdata(ng).q = q;

  xcenter = xlow + dx/2 + (0:(mx-1))*dx;
  xedge = xcenter(1:mx) + dx/2;

  % for compatibility with old matlab41/plotframe1 convention:
  x = xcenter;

  qmin = min([min(q), qmin]);
  qmax = max([max(q), qmax]);

  nplots = size(q,2);
  if (mappedgrid == 1)
    xp = mapc2p(xcenter);
  else
    xp = xcenter;
  end;

  if (nplots > 1)
    for n = 1:nplots,
      subplot(nplots,1,n)
      newplot;
      add_line2plot(xp,q(:,n),level,markerstyle{level},...
	  linecolors{level},linestyle{level});
    end;
  else
    add_line2plot(xp,q,level,markerstyle{level},...
	linecolors{level},linestyle{level});
  end;

if exist('aftergrid')==2
  % make an m-file with this name for any other commands you
  % want executed at the end of drawing each grid
  aftergrid;
end;

end  % loop on ng

% add title and labels:
if UserVariable == 1
  str = sprintf('%s at time %8.4f',UserVariableFile,t);
  title(str,'fontsize',15);
elseif (nplots == 1)
  str = sprintf('q(%d) at time %8.4f',mq,t);
  title(str,'fontsize',15);
else
  for  n = 1:nplots,
    subplot(nplots,1,n);
    str = sprintf('q(%d) at time %8.4f',n,t);
    title(str);
  end  % loop on mq
end  % if UserVariable

if exist('afterframe')==2
  % make an m-file with this name for any other commands you
  % want executed at the end of drawing each frame
  % for example to change the axes, or add something to the plot

  afterframe;
end
