% 3D Graphics.
% PlotType = 4 : Scatter plots for spherically symmetric data
%                Specify (x0,y0,z0) in setplot3.m
%
% PlotType = 5 : Isosurfaces, with lighting.  Specify surface level
%                in setplot3.m
%
% PlotType = 6 : Slice values.  Specify slice coordinates in setplot3.m.
%

% read time and number of grids:
n1 = Frame+10000;
fname = ['fort.',num2str(n1)];
fname(6) = 't';

if ~exist(fname)
  disp(' ')
  disp(['Frame ',num2str(Frame),' does not exist ***']);
  disp(' ')
end

% Check for proper input

% Do some checking to make sure input is right..
if (PlotType <= 3)
  % ---------
  % Slices
  % ---------
  global patch_handles patch_count;
  global contour_handles;
  if (PlotType == 3)
    % Change color map for schlieren plot
    colormap(flipud(gray(2048)).^5)

    if (~isempty(ContourValues))
      fprintf('\n');
      fprintf(' *** Warning : ContourValues is set to a non-empty matrix.\n');
      fprintf('     Schlieren data will be contoured.\n\n');
    end;
  end;

  % In case the slice and contour values just got commented
  % out of setplot3.
  if (~exist('ContourValues'))
    ContourValues = [];
  end;
  if (~exist('xSliceCoords'))
    xSliceCoords = [];
  end;
  if (~exist('ySliceCoords'))
    ySliceCoords = [];
  end;
  if (~exist('zSliceCoords'))
    zSliceCoords = [];
  end;

  xpatch_count = zeros(length(xSliceCoords),MaxLevels);
  ypatch_count = zeros(length(ySliceCoords),MaxLevels);
  zpatch_count = zeros(length(zSliceCoords),MaxLevels);

  % 3rd index used to store multiple patches at same level
  xpatch_handles = zeros(length(xSliceCoords),MaxLevels,1);
  ypatch_handles = zeros(length(ySliceCoords),MaxLevels,1);
  zpatch_handles = zeros(length(zSliceCoords),MaxLevels,1);

  % We use cell arrays here because each patch will have in
  % general a vector of line objects associated with it, if
  % the user has specified contour values.  Cell arrays are
  % easier to use here than a four dimensional array.
  % The third dimension is used to store multiple patches at
  % same level, on same slice.
  xcontour_handles = cell(length(xSliceCoords),MaxLevels,1);
  ycontour_handles = cell(length(ySliceCoords),MaxLevels,1);
  zcontour_handles = cell(length(zSliceCoords),MaxLevels,1);

elseif PlotType == 4
  % -------------
  % Scatter plots
  % -------------
  maxr = 0;

  if (size(PlotStyle,1)) == 1 % scalar type provided
    % expand it so we can use it below.
    sym = PlotStyle(1,:);
    plotstyle = sym;
    for i = 1:MaxLevels-1,
      plotstyle = strvcat(plotstyle,sym);
    end;
  else
    plotstyle = PlotStyle;
  end;

elseif (PlotType == 5)
  % isosurfaces....

  if sum(PlotData > 0) > 1 & ngrids > 1
    fprintf('\n');
    fprintf(' *** Warning: Isosurfaces are being plotted for more than\n');
    fprintf('     one level of AMR data.  You may want to set PlotData so \n');
    fprintf('     only one value is nonzero.\n\n');
  end;

  versionstring = version;
  version_number = str2num(versionstring(1));
  if (version_number < 6)
    fprintf('\n');
    fprintf(' *** Problem in plotframe3: Transparency supported only in\n')
    fprintf('     Matlab versions 6 and higher.\n\n');
  end

  if (sum(IsosurfAlphas) ~= length(IsosurfAlphas))
    % Not all transparencies are set to 1 (=opaque).
    % User has defined some surfaces to be transparent
    % OpenGL Renderer is required.

    rset = set(gcf,'Renderer');
    found_opengl = 0;
    for i = 1:length(rset),
      if (strcmp(rset(i),'OpenGL'))
	found_opengl = 1;
	break;
      end;
    end;

    if (~found_opengl)
      fprintf(' *** Warning : The graphics renderer OpenGL does not \n');
      fprintf('	              appear to be available on your system.\n');
      fprintf('		      Surfaces will not appear transparent.\n\n');
    else
      set(gcf,'Renderer','OpenGL');
    end;
  end;

  % Check length of IsosurfColors and IsosurfAlphas
  if (length(IsosurfAlphas) < length(IsosurfValues))
    error(['Size of IsosurfAlphas must be >= size of IsosurfValues.']);
  elseif (size(IsosurfColors,1) < length(IsosurfValues))
    error(['Size of IsosurfColors must be >= size of IsosurfValues.']);
  end;

end;


% End of input check.

% If there is a file 'fort.tXXXX'
if exist(fname)

  % start reading data, beginning with parameters for this frame:

  % Read data from fname = 'fort.tXXXX'
  fid = fopen(fname);

  t = fscanf(fid,'%g',1);        fscanf(fid,'%s',1);
  meqn = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
  ngrids = fscanf(fid,'%d',1);   fscanf(fid,'%s',1);
  fclose(fid);

  disp(' ');
  disp(['Frame ',num2str(Frame),' at time t = ',num2str(t)]);

  qmin = [];
  qmax = [];
  ncells = [];

  clf;
  if exist('beforeframe')==2
    % make an m-file with this name for any other commands you
    % want executed before drawing each frame, for example
    % if you want to use axes to specify exactly where the
    % plot will be in the window, aspect ratio, etc.
    beforeframe;
  end

  % change the file name to read the q data:
  fname(6) = 'q';
  fid = fopen(fname);


  %=============================================
  % MAIN LOOP ON GRIDS FOR THIS FRAME:
  %=============================================

  for ng = 1:ngrids

    % read parameters for this grid:

    gridno = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    level = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    mx = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    my = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    mz = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);

    xlow = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    ylow = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    zlow = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);

    dx = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    dy = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
    dz = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);

    % Get grid info for plotting grid boxes.
    xhigh = xlow + mx*dx;
    yhigh = ylow + my*dy;
    zhigh = zlow + mz*dz;

    xe = xlow:dx:(xhigh+1e-8);
    ye = ylow:dy:(yhigh+1e-8);
    ze = zlow:dz:(zhigh+1e-8);

    xe(mx+1) = xhigh;
    ye(my+1) = yhigh;
    ze(mz+1) = zhigh;

    xc = xe(1:mx) + dx/2;
    yc = ye(1:my) + dy/2;
    zc = ze(1:mz) + dz/2;


    [xcm,ycm,zcm] = meshgrid(xc,yc,zc);

    if (PlotCubeEdges(level) == 1 & PlotType <= 3)
      plotCube(xlow,xhigh,ylow,yhigh,zlow,zhigh);
    end;

    % read q data:

    data = fscanf(fid,'%g',[meqn,mx*my*mz]);

    % if we're not plotting data at this level, skip to next grid
    if PlotData(level) == 1

      data = data';

      if (UserVariable == 1)
	% User has supplied a function to convert original q variables to
	% the variable which is to be plotted, e.g. Mach number, entropy.
	qdata = feval(UserVariableFile,xc,yc,zc,data);
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
	% Slices

        if length(ContourValues)==1
           disp(['** Warning: Contour values will be chosen independently'...
                 ' on each slice'])
           disp('   Set ContourValues to a vector of values to insure they match')
           end

	% The routine intersectQWithSlices uses the global data
	% patch_handles and patch_count and changes their values.

	% Now plot any x slices.
	if (~isempty(xSliceCoords))
	  sdir = 'x';
	  patch_handles = xpatch_handles;
	  patch_count = xpatch_count;
	  contour_handles = xcontour_handles;
	  intersectQWithSlices(xc,yc,zc,xe,ye,ze,qmesh,sdir,...
	      xSliceCoords,level,PlotGrid(level),...
	      PlotGridEdges(level),PlotType,ContourValues);
	  xcontour_handles = contour_handles;
	  xpatch_handles = patch_handles;
	  xpatch_count = patch_count;
	end;

	if (~isempty(ySliceCoords))
	  sdir = 'y';
	  patch_handles = ypatch_handles;
	  patch_count = ypatch_count;
	  contour_handles = ycontour_handles;
	  intersectQWithSlices(xc,yc,zc,xe,ye,ze,qmesh,sdir,...
	      ySliceCoords,level,PlotGrid(level),...
	      PlotGridEdges(level),PlotType,ContourValues);
	  ycontour_handles = contour_handles;
	  ypatch_handles = patch_handles;
	  ypatch_count = patch_count;
	end;

	if (~isempty(zSliceCoords))
	  sdir = 'z';
	  patch_handles = zpatch_handles;
	  patch_count = zpatch_count;
	  contour_handles = zcontour_handles;
	  intersectQWithSlices(xc,yc,zc,xe,ye,ze,qmesh,sdir,...
	      zSliceCoords,level, PlotGrid(level),...
	      PlotGridEdges(level),PlotType,ContourValues);
	  zcontour_handles = contour_handles;
	  zpatch_handles = patch_handles;
	  zpatch_count = patch_count;
	end;

      elseif PlotType == 4
	% Scatter plot for radially symmetric data.

	r = sqrt((xcm - x0).^2 + (ycm - y0).^2 + (zcm - z0).^2);
	rvec = reshape(r,prod(size(r)),1);
	qvec = reshape(qmesh,prod(size(qmesh)),1);
	% note that plotstyle has been expanded from above, if it was a
	% scalar value.
	plot(rvec,qvec,plotstyle(level,:));
	maxr = max([maxr; rvec]);
	view(2);
	% set(gca,'XLim',[0 maxr]);
	hold on;

      elseif PlotType == 5
	% Isosurfaces

	% Loop over surface values.
	for i = 1:length(IsosurfValues),
	  hsurf = patch(isosurface(xcm,ycm,zcm,qmesh,IsosurfValues(i)));
	  isonormals(xcm,ycm,zcm,qmesh,hsurf);
	  set(hsurf,'FaceColor',IsosurfColors(i,:),'EdgeColor','none');
	  if (version_number >= 6)
	    % Property 'FaceAlpha' supported only in versions 6 and higher.
	    set(hsurf,'FaceAlpha',IsosurfAlphas(i));
	  end;
	  daspect([1 1 1]);
	  view(3);
	  camlight;
	  lighting phong;
	  grid on;
	  hold on;
	end

      end; % end PlotType Choice

      % -----------------------------------------------

    end %   if PlotData(level)==1
  %=============================================
  end % loop on ng (plot commands for each grid)
  %=============================================

  % add title and labels:

  if UserVariable == 1
	title([UserVariableFile,'   at time t = ', num2str(t)])
      else
	title(['q(',num2str(mq),')   at time t = ', num2str(t)])
      end
 
 if (PlotType ~= 4)
   xlabel('x')
   ylabel('y')
   zlabel('z')
   end


 if (PlotType <= 3)
   % Now plot intersections of x-y planes, x-z planes and y-z planes.
   v2 = ones(1,2);

   xyIntersect = cell(length(xSliceCoords),length(ySliceCoords));
   xzIntersect = cell(length(xSliceCoords),length(zSliceCoords));
   yzIntersect = cell(length(ySliceCoords),length(zSliceCoords));

   for level = 1:length(PlotData),
     if PlotData(level) == 1
       if (~isempty(xSliceCoords) & ~isempty(ySliceCoords))
	 for i = 1:length(xSliceCoords),
	   for j = 1:length(ySliceCoords),
	     for kx = 1:xpatch_count(i,level)
	       for ky = 1:ypatch_count(j,level)
		 udatax = get(xpatch_handles(i,level,kx),'UserData');
		 udatay = get(ypatch_handles(j,level,ky),'UserData');
		 % Two patch slices might intersect
		 xval = udatax.sval;
		 yval = udatay.sval;
		 zmin = udatax.zmin; % Should be the same as udatay.zmin
		 zmax = udatax.zmax; %     "       "         udatay.zmax
		 if (udatax.ymin > yval | udatax.ymax < yval)
		   % no intersection
		 elseif (udatay.xmin > xval | udatay.xmax < xval)
		   % also no intersection
		 else
		   h = line('XData',xval*v2,'YData',yval*v2,...
		       'ZData',[zmin zmax],'Color','k');
		   set(h,'Tag','xyintersect');
		   xyIntersect{i,j} = [xyIntersect{i,j}, h];
		 end; % intersection check
	       end; % ky loop
	     end; % kx loop
	   end; % yslice(j) loop
	 end; % xslice(i) loop
       end; % isempty loop

       if (~isempty(xSliceCoords) & ~isempty(zSliceCoords))
	 for i = 1:length(xSliceCoords),
	   for j = 1:length(zSliceCoords),
	     for kx = 1:xpatch_count(i,level)
	       for kz = 1:zpatch_count(j,level)
		 udatax = get(xpatch_handles(i,level,kx),'UserData');
		 udataz = get(zpatch_handles(j,level,kz),'UserData');
		 % Two patch slices might intersect
		 xval = udatax.sval;
		 zval = udataz.sval;
		 ymin = udatax.ymin; % Should be the same as udatay.zmin
		 ymax = udatax.ymax; %     "       "         udatay.zmax
		 if (udatax.zmin > zval | udatax.zmax < zval)
		   % no intersection
		 elseif (udataz.xmin > xval | udataz.xmax < xval)
		   % also no intersection
		 else
		   h = line('XData',xval*v2,'YData',[ymin ymax],...
		       'ZData',zval*v2,'Color','k');
		   set(h,'Tag','xzintersect');
		   xzIntersect{i,j} = [xzIntersect{i,j}, h];
		 end; % intersection check
	       end; % kz loop
	     end; % kx loop
	   end; % zslice(j) loop
	 end; % xslice(i) loop
       end; % isempty loop

       if (~isempty(ySliceCoords) & ~isempty(zSliceCoords))
	 for i = 1:length(ySliceCoords),
	   for j = 1:length(zSliceCoords),
	     for ky = 1:ypatch_count(i,level)
	       for kz = 1:zpatch_count(j,level)
		 udatay = get(ypatch_handles(i,level,ky),'UserData');
		 udataz = get(zpatch_handles(j,level,kz),'UserData');
		 % Two patch slices might intersect
		 yval = udatay.sval;
		 zval = udataz.sval;
		 xmin = udatay.xmin; % Should be the same as udatay.zmin
		 xmax = udatay.xmax; %     "       "         udatay.zmax
		 if (udatay.zmin > zval | udatay.zmax < zval)
		   % no intersection
		 elseif (udataz.ymin > yval | udataz.ymax < yval)
		   % also no intersection
		 else
		   h = line('XData',[xmin xmax],'YData',yval*v2,...
		       'ZData',zval*v2,'Color','k');
		   set(h,'Tag','yzintersect');
		   yzIntersect{i,j} = [yzIntersect{i,j}, h];
		 end; % intersection check
	       end; % kz loop
	     end; % ky loop
	   end; % zslice(j) loop
	 end; % yslice(i) loop
       end; % isempty loop

     end; % We plotted data at this level
   end; % end level loop
 end;

 % done with this frame
 status = fclose(fid);

 if exist('afterframe')==2
   afterframe  % make an m-file with this name for any other commands you
   % want executed at the end of drawing each frame
   % for example to change the axes, or add a curve for a
   % boundary
 end;

 hold off;

end % if exist(fname)
