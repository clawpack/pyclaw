%
% PLOTFRAME2 plots data from a single Clawpack output file.
%
%    PLOTFRAME2 is called from PLOTCLAW2, the driver script for Clawpack
%    graphics.  PLOTFRAME2 uses parameters, defined in the user's
%    workspace, usually in the SETPLOT2 script, to determine what kind of
%    plot to create, and various other plotting options.
%
%    See SETPLOT for a complete list of the parameters that the user can
%    set.
%
%    There are two basic types of plots available.  These are
%
%         Surface plots : 2d plots of the data.  This data may be viewed as
%         a manifold, if desired.
%
%         Scatter plots : Spherically symmetric data is plotted as a
%         function of a 1d variable, usually the distance from a fixed point
%         (x0,y0).  The result is a line-type plot of points representing the
%         data.
%
%    See SETPLOT, PLOTCLAW2, MANIFOLD.

if (~exist('amrdata'))
  error('*** plotframe2 : ''amrdata'' does not exist;  Call readamrdata or plotclaw2');
end;

if length(amrdata) == 0
  disp(' ')
  disp(['Frame ',num2str(Frame),'(',outputflag,') does not exist ***']);
  disp(' ')
  return;
end

disp(' ');
disp(['Frame ',num2str(Frame),' at time t = ',num2str(t)]);

if exist('beforeframe')==2
  % make an m-file with this name for any other commands you
  % want executed before drawing each frame, for example
  % if you want to use axes to specify exactly where the
  % plot will be in the window, aspect ratio, etc.
  beforeframe;
end


% Do some checking to make sure input is right..
if (PlotType <= 3)

  set_value('mappedgrid','MappedGrid',0);
  if (mappedgrid == 1 & ~exist('mapc2p'))
    error('*** MappedGrid = 1 but no ''mapc2p'' function was found.');
  end;

  set_value('manifold','Manifold',0);
  if (manifold == 1 & ~exist('mapc2m'))
    error('*** Manifold = 1, but no ''mapc2m'' function was found.');
  end;

  if (manifold == 1)
    set_value('view_arg','UserView',3);
  else
    set_value('view_arg','UserView',2);
  end;

  set_value('cvalues','ContourValues',[]);
  if (length(cvalues) == 1 & length(amrdata) > 1)
    disp(' ');
    disp(' *** Warning: Contour values will be chosen independently');
    disp(' *** on each mesh patch. Set ContourValues to a vector of values');
    disp(' *** to insure contour lines match across mesh boundaries.');
  end;


  if (PlotType == 2)
    if (isempty(cvalues))
      disp(' ');
      disp(' *** Warning : You have specified PlotType == 2, but have');
      disp(' *** not set any contour values. Set ContourValues to a ');
      disp(' *** vector of values or integer (number of contour lines). ');
    end;
  end;


  if (PlotType == 3)
    % Change color map for schlieren plot
    % colormap(flipud(gray(2048)).^5)
    colormap(flipud(gray(2048)).^10)

    if (~isempty(cvalues))
      disp(' ');
      disp(' *** Warning : ContourValues is set to a non-empty matrix.');
      disp(' *** Schlieren data will be contoured.');
    end;
  end;

  xscoords = [];
  yscoords = [];
  zscoords = 0;

  sliceCoords = {xscoords, yscoords, zscoords};
  sdir = {'x','y','z'};

  % Stored in UserData of current figure (i.e. gcf);
  clear_amrplot;
  create_amrplot(MaxLevels,xscoords, yscoords, zscoords);
  newplot;  % Use this to deal with any hold on/ hold off issues.

end; % end of input checking for PlotType <= 3

if PlotType == 4
  % -------------
  % Scatter plots
  % -------------

  view_arg = 2;

  set_value('usermap1d','UserMap1d',0);
  if (usermap1d == 1)
    if (~exist('map1d'))
      error('*** You have set UserMap1d=1, but no ''map1d'' function was found');
    end;
  else
    if ((~exist('x0') | ~exist('y0')))
      str = sprintf(['*** plotframe2 : (x0,y0) must be defined before you\n',...
	    '*** can do a line or scatter plot, or set UserMap1d and create\n',...
	    '*** a function ''map1d''.']);
      error(str);
    end;
  end;

  set_value('mappedgrid','MappedGrid',0);

  if (exist('ScatterStyle'))
    pstyle = ScatterStyle;
  elseif (exist('LineStyle'))
    pstyle = LineStyle;
  else
    error([' *** plotframe2 : Set either ''ScatterStyle'' or ',...
	  '''LineStyle''.']);
  end;

  if (~iscell(pstyle))
    error(['*** plotframe2 : ScatterStyle or LineStyle must be',...
	  'cell matrices.  Use ''setplotstyle'' to set either of these ',...
	  'variables']);
  end;

  [linestyle,linecolors,markerstyle] = get_plotstyle(pstyle,MaxLevels);

  clear_amrplot;
  create_amrplot(MaxLevels);
  newplot;   % deals with some hold on/off issues...

end;  % end of input checking for PlotType == 4

qmin = [];
qmax = [];
ncells = [];

%=============================================
% MAIN LOOP ON GRIDS FOR THIS FRAME:
%=============================================

ngrids = length(amrdata);  % length of structure array
for ng = 1:ngrids,

  gridno = amrdata(ng).gridno;
  level = amrdata(ng).level;

  % if we're not plotting data at this level, skip to next grid
  if (PlotData(level) == 0)
    continue;
  end;

  % Set block number for multi-block calculations.
  set_blocknumber(gridno);

  mx = amrdata(ng).mx;
  my = amrdata(ng).my;

  xlow = amrdata(ng).xlow;
  ylow = amrdata(ng).ylow;

  dx = amrdata(ng).dx;
  dy = amrdata(ng).dy;

  xedge = xlow + (0:mx)*dx;
  yedge = ylow + (0:my)*dy;

  xcenter = xedge(1:mx) + dx/2;
  ycenter = yedge(1:my) + dy/2;

  % for compatibility with old matlab41/plotframe2 convention:
  x = xcenter;
  y = ycenter;

  % read q data:
  data = amrdata(ng).data;


  data = data';

  if (UserVariable == 1)
    % User has supplied a function to convert original q variables to
    % the variable which is to be plotted, e.g. Mach number, entropy.
    qdata = feval(UserVariableFile,data);
    q = reshape(qdata,mx,my);
  else
    q = reshape(data(:,mq),mx,my);
  end

  amrdata(ng).q = q;

  % q must be permuted so that it matches the dimensions of xcm,ycm,zcm
  % created by meshgrid.
  if (PlotType == 3)
    % Scheieren plot;  we plot the gradient, not the values.
    [qx,qy] = gradient(q,dx,dy);
    qs = sqrt(qx.^2 + qy.^2);
    qmesh = qs';
  else
    qmesh = q';
  end;

  % minimum over all grids at this time, but not necessarily on slice
  % shown.
  qmin = min([qmin,min(min(q))]);
  qmax = max([qmax,max(max(q))]);

  % keep count of how many cells at this refinement level:
  if length(ncells) < level
    ncells(level) = 0;
  end;
  ncells(level) = ncells(level) + mx*my;

  % -----------------------------------------------
  % plot commands go here
  if PlotType <= 3

    % Add amr patch of manifold into current plot.
    sval = level - 6;  % for top down viewing
    zedge = [sval sval];
    zcenter = [sval sval];
    sdir = 'z';
    snum = 1;   % only one slice in 2d plot
    % only mask patches underneath if we are plotting a Manifold
    maskflag = (manifold == 1);
    add_patch2slice(sdir,sval,snum,xcenter,ycenter,zcenter, ...
	xedge,yedge,zedge,qmesh,level,...
	cvalues,mappedgrid,manifold,maskflag,ng);
  end;  % end of plotting for PlotType == 3

  if (PlotType == 4)
    % 1d Line plots

    [xcm,ycm] = meshgrid(xcenter,ycenter);

    if (usermap1d == 1)
      if (mappedgrid == 1)
	[xpm,ypm] = mapc2p(xcm,ycm);
        [rvec,qvec] = map1d(xpm,ypm,qmesh);
      else
        [rvec,qvec] = map1d(xcm,ycm,qmesh);
      end
      [rs,cs] = size(rvec);
      [rq,cq] = size(qvec);
      if (cs > 1 | cq > 1)
	error(['plotframe2 : map1d can only return single columns vectors ',...
	       'for s or q']);
      end;
    else
      if (mappedgrid == 1)
	[xpm,ypm] = mapc2p(xcm,ycm);
	r = sqrt((xpm - x0).^2 + (ypm - y0).^2);
      else
	r = sqrt((xcm-x0).^2 + (ycm - y0).^2);
      end
      rvec = reshape(r,prod(size(r)),1);
      qvec = reshape(qmesh,prod(size(qmesh)),1);
    end;

   add_line2plot(rvec,qvec,level,markerstyle{level},...
	linecolors{level},linestyle{level});

  end; % end of plotting for PlotType == 4

if exist('aftergrid')==2
  % make an m-file with this name for any other commands you
  % want executed at the end of drawing each grid
  aftergrid;
end;

end % loop on ng (plot commands for each grid)
%=============================================

% Set user-defined variables from setplot2.m :

% add title and labels:
if UserVariable == 1
  str = sprintf('%s at time %8.4f',UserVariableFile,t);
  title(str,'fontsize',15);
else
  str = sprintf('q(%d) at time %8.4f',mq,t);
  title(str,'fontsize',15);
end

if (PlotType <= 3)
  setPlotGrid(PlotGrid);
  setPlotGridEdges(PlotGridEdges);
  if (PlotType == 2)
    setslicecolor('w');
  end;
  if (manifold == 0 & view_arg == 2)
    % To really make sure we get a 2d view and can see all levels.
    set(gca,'ZLimMode','auto');
  end;
  xlabel('x')
  ylabel('y')
end;

% Set view point
view(view_arg);

if exist('afterframe')==2
  % make an m-file with this name for any other commands you
  % want executed at the end of drawing each frame
  % for example to change the axes, or add a curve for a
  % boundary

  afterframe;
end;
