
% claw/matlab/plotframe1.m
%
% Plot data at a single time produced by clawpack, 
%    from fort.t**** and fort.q**** files at this time frame
% Assumes Frame is already set.  
% Usually called from plotclaw1.m within a loop on Frame, but
% can be used alone to plot one specific Frame.
%
% Assumes variables described in claw/matlab/setplot1.m are set, e.g. 
% by calling that routine.

% Modified 10/01 for bearclaw to allow multiple grids in each frame.
% The data xgrid and datagrid from each grid is read in and combined into 
% single arrays xmulti and datamulti that contains the finest level
% information available at each point.  Data manipulation and plot commands 
% are then applied to these arrays.

% The data from each grid separately is also stored into arrays
% xgridN and datagridN where N has the value 1..ngrids.   These are not used
% directly but are available if you want to examine all the data.
% For example, the commands
%   plot(xgrid2,datagrid2(:,1),'o')
%   hold on
%   plot(xgrid3,datagrid3(:,1),'+')
%   hold off
% would plot all data values q(1) from grids 2 and 3 together.

% Type i at the query prompt to get information about each grid, its
% level, xlower, mx, and dx.


% read time and number of grids:

n1 = Frame+10000;
fname = ['fort.',num2str(n1)];
fname(6) = 't';
 
if ~exist(fname) 
     disp(' ');
     disp(['Frame ',num2str(Frame),' does not exist ***']);
     end

if exist(fname)

 % start reading data, beginning with parameters for this frame:

 fid = fopen(fname);

 t = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
 meqn = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
 ngrids = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);

 if ngrids==0 
   % flag that the computation has ended, in case nframes wasn't known
   % initially
   break
   end

 disp(' ')
 disp(['Frame ',num2str(Frame),' at time t = ',num2str(t)]);


 if exist('beforeframe')==2
    beforeframe  % make an m-file with this name for any other commands you
                % want executed before drawing each frame, for example
                 % if you want to use axes to specify exactly where the
                 % plot will be in the window, aspect ratio, etc.
    end


 % change the file name to read the q data:
 fname(6) = 'q';
 fid = fopen(fname);

 % initialize arrays in which to accumulate full domain solution
 % by piecing together data from each grid:
 xmulti = [];
 datamulti = [];

 %=============================================
 % MAIN LOOP ON GRIDS FOR THIS FRAME:
 %=============================================

 for ng = 1:ngrids

   % read parameters for this grid:

   gridno = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
   level = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
   mx = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);

   xlow = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
   dx = fscanf(fid,'%g',1);     fscanf(fid,'%s',1);

   GridInfo(ng,1) = gridno;
   GridInfo(ng,2) = level;
   GridInfo(ng,3) = mx;
   GridInfo(ng,4) = xlow;
   GridInfo(ng,5) = dx;

   % read q data:

   data = fscanf(fid,'%g',[meqn,mx]);
   datagrid = data';
   xgrid = xlow + ((1:mx) - 0.5)*dx;
   dataname = ['datagrid' num2str(ng)];
   xname = ['xgrid' num2str(ng)];
   eval([dataname '= datagrid;']);
   eval([xname '= xgrid;']);

   % insert this data in place of any coarser grid data already in 
   % xmulti and datamulti:

   xup = xlow + mx*dx;
   ileft = find(xmulti < xlow);  % need to keep coarse grid data here
   iright = find(xmulti > xup);  % need to keep coarse grid data here
   xmulti = [xmulti(ileft)  xgrid  xmulti(iright)];
   datamulti = [datamulti(ileft,:);
                datagrid;
                datamulti(iright,:)]; 
   end  % loop on ng

   % rename full domain data for use in usual plot commands below:
   x = xmulti;
   data = datamulti;

   if UserVariable==1
       % User has supplied a function to convert original q variables to
       % the variable which is to be plotted, e.g. Mach number, entropy.
       q = feval(UserVariableFile,data);
       end


   if MappedGrid
      % coordinate mapping must be applied
      x = mapc2p(x);
      end

 %=====================================================================
 %  The plot command:
 %=====================================================================

 if UserVariable==1
     % q has been computed using user's file.
     plot(x,q,plotstyle)
     title([UserVariableFile,'   at time t = ', num2str(t)])
  elseif length(mq)==1
     % plotting a single variable
     q = data(:,mq);
     plot(x,q,plotstyle)
     title(['q(',num2str(mq),')   at time t = ', num2str(t)])
  else
     % plotting more than one variable in same figure
     % put each in a separate subplot
     for  mq1=mq
        subplot(length(mq),1,find(mq==mq1))
        q = data(:,mq1);
        plot(x,q,plotstyle)
        title(['q(',num2str(mq1),')   at time t = ', num2str(t)])
        end  % loop on mq1
     end  % if UserVariable

 status = fclose(fid);

 if exist('afterframe')==2  
    afterframe  % make an m-file with this name for any other commands you
	        % want executed at the end of drawing each frame
                % for example to change the axes, or add something to the plot
    end


end % if exist(fname)
