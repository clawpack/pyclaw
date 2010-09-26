%
% PLOTFRAME3 plots data from a single Clawpack output file.
%
%    PLOTFRAME3 is called from PLOTCLAW3, the driver script for Clawpack
%    graphics.  PLOTFRAME3 uses parameters, defined in the user's
%    workspace, usually in the SETPLOT3 script, to determine what kind of
%    plot to create, and various other plotting options.
%
%    See SETPLOT for a complete list of the parameters that the user can
%    set.
%
%    There are three basic types of plots available.  These are
%
%         Slice plots : 3d slices, along Cartesian or mapped Cartesian
%         coordinates shown.
%
%         Scatter plots : Spherically symmetric data is plotted as a
%         function of a 1d variable, usually the distance from a fixed point
%         (x0,y0,z0).  The result is a 2d line-type plot of points
%         representing the data.
%
%         Isosurfaces : Surfaces of constant data values are shown. These
%         can be combined with slices, if desired.
%
%    See SETPLOT, PLOTCLAW3, MAPPEDGRID.

if (~exist('amrdata'))
  str = sprintf(['*** Plotframe3 : ''amrdata'' does not exist. call\n',...
	         '*** readamrdata or plotclaw3']);
  error(str);
end;

if length(amrdata) == 0
  disp(' ')
  disp(['Frame ',num2str(Frame),'(',outputflag,') does not exist ***']);
  disp(' ')
  return;
end;

disp(' ');
disp(['Frame ',num2str(Frame),' at time t = ',num2str(t)]);

if exist('beforeframe')==2
  % make an m-file with this name for any other commands you
  % want executed before drawing each frame, for example
  % if you want to use axes to specify exactly where the
  % plot will be in the window, aspect ratio, etc.
  beforeframe;
end

% ------------------------------------------------------------------------
% Check PlotType <= 3 inputs (slices)
% ------------------------------------------------------------------------
if (PlotType <= 3)

  % The manifold option isn't used in 3d
  manifold = 0;

  set_value('mappedgrid','MappedGrid',0);
  set_value('view_arg','UserView',3);
  set_value('cvalues','ContourValues',[]);
  set_value('interpmethod','InterpMethod','*linear');

  if (mappedgrid == 1)
    if (~exist('mapc2p'))
      error('*** MappedGrid == 1 but no ''mapc2p'' function was found');
    end;
  end;

  if (PlotType == 2)
    if (isempty(cvalues))
      disp(' ');
      disp(' *** Warning : You have specified PlotType == 2, but have');
      disp(' *** not set any contour values. Set ContourValues to a ');
      disp(' *** vector of values or integer (number of contour lines). ');
    end;
  end;

  % Change color map for schlieren plot
  if (PlotType == 3)

    % colormap(flipud(gray(2048)).^5)
    colormap(flipud(gray(2048)).^10)

    if (~isempty(cvalues))
      disp(' ');
      disp(' *** Warning : ContourValues is set to a non-empty matrix.');
      disp(' *** Schlieren data will be contoured.');
    end;
  end;

  set_value('xscoords','xSliceCoords',[]);
  set_value('yscoords','ySliceCoords',[]);
  set_value('zscoords','zSliceCoords',[]);

  if (length(cvalues) == 1 & (length(amrdata) > 1 | ...
    (length(xscoords) + length(yscoords) + length(zscoords) > 1)))
    disp(' ');
    disp(' *** Warning: Contour values will be chosen independently');
    disp(' *** on each amr patch or slice. Set ContourValues to a ');
    disp(' *** vector of values to insure that contour lines match ');
    disp(' *** across mesh boundaries or slices.');
  end;

  set_value('isosurfvalues','IsosurfValues',[]);

  m = length(isosurfvalues);
  if (m > 0)
    % now check input for isosurfaces

    set_value('isoalphas','IsosurfAlphas',1);
    isosurfalphas = set_length(isoalphas,m);

    set_value('isc','IsosurfColors','q');
    if (iscell(isc))
      % In case SETPLOTSTYLE was used to to set the isosurfcolors.
      isocolors = strvcat(isc);
    else
      isocolors = isc;
    end;
    isosurfcolors = set_length(isocolors,m);

    versionstring = version;
    version_number = str2num(versionstring(1));

    for i = 1:length(isosurfalphas)
      a = isosurfalphas(i);
      has_opengl = strcmp(get(gcf,'Renderer'),'OpenGL');
      if (a < 1 & version_number < 6)
	disp(' ');
	disp(' *** plotframe3: Transparency for isosurfaces supported only in')
	disp(' *** Matlab versions 6 and higher.');
	isosurfalphas(i) = 1;
      elseif (a < 1 & has_opengl == 0)
	disp(' *** Warning : The graphics renderer OpenGL is not set on your');
	disp(' *** system.  Use setopengl to set this renderer.  Without it,');
	disp(' *** isosurfaces will not appear transparent');
	isosurfalphas(i) = 1;
      end;
    end;
  else
    % These won't be used...
    isosurfcolors = '';
    isosurfalphas = [];
  end;  % end of length(isosurfvalues > 0)

  clear_amrplot;
  create_amrplot(MaxLevels,xscoords,yscoords,zscoords,isosurfvalues);
  newplot;

end;  % end of PlotType <= 3

% ------------------------------------------------------------------------
% Check PlotType == 4 inputs (1d line plots)
% ------------------------------------------------------------------------
if (PlotType == 4)

  view_arg = 2;

  set_value('usermap1d','UserMap1d',0);
  if (usermap1d == 1)
    if (~exist('map1d'))
      error('*** You have set UserMap1d=1, but no ''map1d'' function was found');
    end;
  else
    if ((~exist('x0') | ~exist('y0') | ~exist('z0')))
      str = sprintf(['*** plotframe3 : (x0,y0,z0) must be defined before you\n',...
	    '*** can do a line or scatter plot.  Or, create a ''map1d''',...
	    ' function and set UserMap1d = 1']);
      error(str);
    end;
  end;

  set_value('mappedgrid','MappedGrid',0);

  if (exist('ScatterStyle'))
    pstyle = ScatterStyle;
  elseif (exist('LineStyle'))
    pstyle = LineStyle;
  else
    error([' *** plotframe3 : Set either ''ScatterStyle'' or ',...
	  '''LineStyle''.']);
  end;

  if (~iscell(pstyle))
    str = sprintf(['*** plotframe3 : ScatterStyle or LineStyle must be cell \n',...
	  'matrices.  Use SETPLOTSTYLE to set either of these variables.']);
    error(str);
  end;

  [linestyle,linecolors,markerstyle] = get_plotstyle(pstyle,MaxLevels);

  clear_amrplot;
  create_amrplot(MaxLevels);
  newplot;

end; % end of PlotType == 4
% ------------------------------------------------------------------------
% Check PlotType == 5 inputs (deprecated isosurface plots)
% ------------------------------------------------------------------------
if (PlotType == 5)
  disp(' ')
  disp(' *** PlotType = 5 is no longer supported as of Clawpack 4.2')
  disp(' *** Choose different setting for PlotType in setplot3.m')
  disp(' *** Isosurface plots are avaiable with PlotType = 1, 2 or 3')
  disp(' *** by setting IsosurfValues and IsosurfColors')
  disp(' ')
end; % end of PlotType == 5
% End of input check.

% ------------------------------------------------------------------
% Store current handle data structure in the UserData property of the
% current figure, gcf.

% ---------------------------------------------------------------------


% ------------------------------------------------------------------------
% Main loop for this frame.
% ------------------------------------------------------------------------

qmin = [];
qmax = [];
ncells = [];

ngrids = length(amrdata);  % length of structure array
for ng = 1:ngrids,

  gridno = amrdata(ng).gridno;
  level = amrdata(ng).level;

  % bn_plotframeN is a global variable used in mapc2p.m if the grid
  % mapping varies across the blocks.
  set_blocknumber(gridno);

  if (PlotData(level) == 0)
    % Continue to next patch
    continue;
  end;

  mx = amrdata(ng).mx;
  my = amrdata(ng).my;
  mz = amrdata(ng).mz;

  xlow = amrdata(ng).xlow;
  ylow = amrdata(ng).ylow;
  zlow = amrdata(ng).zlow;

  dx = amrdata(ng).dx;
  dy = amrdata(ng).dy;
  dz = amrdata(ng).dz;

  data = amrdata(ng).data;

  % Get grid info for plotting grid boxes.
  xhigh = xlow + mx*dx;
  yhigh = ylow + my*dy;
  zhigh = zlow + mz*dz;

  xedge = xlow + (0:mx)*dx;
  yedge = ylow + (0:my)*dy;
  zedge = zlow + (0:mz)*dz;

  xcenter = xedge(1:mx) + dx/2;
  ycenter = yedge(1:my) + dy/2;
  zcenter = zedge(1:mz) + dz/2;

  % for compatibility with old matlab41/plotframe3 convention:
  x = xcenter;
  y = ycenter;
  z = zcenter;

  data = data';

  if (UserVariable == 1)
    % User has supplied a function to convert original q variables to
    % the variable which is to be plotted, e.g. Mach number, entropy.
    qdata = feval(UserVariableFile,data);
    q = reshape(qdata,mx,my,mz);
  else
    q = reshape(data(:,mq),mx,my,mz);
  end


  % q must be permuted so that it matches the dimensions of xcm,ycm,zcm
  % created by meshgrid.
  if (PlotType == 3)
    % Scheieren plot;  we plot the gradient, not the values.
    [qx,qy,qz] = gradient(q,dx,dy,dz);
    qs = sqrt(qx.^2 + qy.^2 + qz.^2);
    qmesh = permute(qs,[2,1,3]);
  else
    qmesh = permute(q,[2,1,3]);
  end;
  amrdata(ng).q = qmesh;

  % minimum over all grids at this time, but not necessarily on slice
  % shown.
  qmin = min([qmin,min(min(min(q)))]);
  qmax = max([qmax,max(max(max(q)))]);

  % keep count of how many cells at this refinement level:
  if length(ncells) < level
    ncells(level) = 0;
  end;
  ncells(level) = ncells(level) + mx*my*mz;

  % -----------------------------------------------
  % plot commands go here
  if PlotType <= 3

    sliceCoords = {xscoords, yscoords, zscoords};
    sdirs = {'x','y','z'};
    maskflag = 1;  % Mask out patch areas

    for idir = 1:3,  % Loop over directions
      slicevals = sliceCoords{idir}; % user specified slice constants
      sdir = sdirs{idir};

      % Loop over all slices that cube of data qmesh might intersect
      for n = 1:length(slicevals),
	sval = slicevals(n);
	[isect,qcm2] = interp_data_3d(xcenter,ycenter,zcenter,...
	    xedge,yedge,zedge,qmesh,sdir,sval,interpmethod);
	if (isect == 1)  % qmesh intersects a user specified slice.
	  % This command adds a patch to the appropriate slice
	  add_patch2slice(sdir,sval, n,xcenter,ycenter,zcenter,...
	      xedge,yedge,zedge,qcm2,level,...
	      cvalues,mappedgrid, manifold,maskflag,ng);
	end;
      end;
    end;

    % Add cube to plot - whether it is visible or not depends on how user
    % set PlotCubeEdges.
    add_cube2plot(xedge,yedge,zedge,level,mappedgrid);

    % Add isosurfaces to plot
    add_isosurface2plot(xcenter,ycenter,zcenter,q,level,isosurfvalues,...
	isosurfcolors,isosurfalphas,mappedgrid);

  end;  % End plotting for PlotType == 3

  if (PlotType == 4)
    % 1d Line plots

    [xcm,ycm,zcm] = meshgrid(xcenter,ycenter,zcenter);

    if (usermap1d == 1)
    if (usermap1d == 1)
      if (mappedgrid == 1)
	[xpm,ypm,zpm] = mapc2p(xcm,ycm,zcm);
        [rvec,qvec] = map1d(xpm,ypm,zpm,qmesh);
      else
        [rvec,qvec] = map1d(xcm,ycm,zcm,qmesh);
      end
      [rs,cs] = size(rvec);
      [rq,cq] = size(qvec);
      if (cs > 1 | cq > 1)
	error(['plotframe3 : map1d can only return single columns vectors ',...
	       'for s or q']);
      end;
    else
      if (mappedgrid == 1)
	[xpm,ypm,zpm] = mapc2p(xcm,ycm,zcm);
	r = sqrt((xpm - x0).^2 + (ypm - y0).^2 + (zpm - z0).^2);
      else
	r = sqrt((xcm - x0).^2 + (ycm - y0).^2 + (zcm - z0).^2);
      end;
      rvec = reshape(r,prod(size(r)),1);
      qvec = reshape(qmesh,prod(size(qmesh)),1);
    end;

    add_line2plot(rvec,qvec,level,markerstyle{level},...
	linecolors{level},linestyle{level});

  end;  % end of plotting for PlotType == 4

if exist('aftergrid')==2
  % make an m-file with this name for any other commands you
  % want executed at the end of drawing each grid
  aftergrid;
end;

end % loop on ng (plot commands for each grid)
% -----------------------------------------------------

set(gca,'FontName','helvetica');
if UserVariable == 1
  str = sprintf('%s at time %8.4f',UserVariableFile,t);
  title(str,'fontsize',15);
else
  str = sprintf('q(%d) at time %8.4f',mq,t);
  title(str,'fontsize',15);
end

% Set user defined grid details, from setplot3.m
if (PlotType <= 3)
  setPlotGrid(PlotGrid);
  setPlotGridEdges(PlotGridEdges);
  setPlotCubeEdges(PlotCubeEdges);
  if (PlotType == 2)
    setslicecolor('w');
  end;

  % Now plot intersections of x-y planes, x-z planes and y-z planes.
  create_sliceintersections(mappedgrid);
  xlabel('x')
  ylabel('y')
  zlabel('z')
end;

view(view_arg);

if exist('afterframe')==2
  afterframe  % make an m-file with this name for any other commands you
  % want executed at the end of drawing each frame
  % for example to change the axes, or add a curve for a
  % boundary
end;
end
