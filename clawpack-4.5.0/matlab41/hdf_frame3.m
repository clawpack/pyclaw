
% first pass at matlab 3d for graphics... plots isosurface.
% this needs to be improved  

 % read time and number of grids:
 n1 = Frame+10000;
 fname = ['q',num2str(n1),'.hdf'];
 fname(2) = '0';
 
 if ~exist(fname) 
     disp(' ');
     disp(['Frame ',num2str(Frame),' does not exist ***']);
%     if Frame==0
%	% no initial data to plot, go on to Frame 1:
%	Frame = 1;
%	fname(length(fname))='1';
%	end
     end

if exist(fname)


% Get some data by reading first grid info 
 SD_id=hdfsd('start',fname,'read');
 [nsds,nattr,status]=hdfsd('fileinfo',SD_id); 
 sds_index=0;
 sds_id=hdfsd('select',SD_id,sds_index);
 [data,status]=hdfsd('readdata',sds_id,0,1,17); 
 status = hdfsd('endaccess',sds_id);
 status = hdfsd('end',SD_id);
 t = data(3);
 meqn = data(4);
 ngrids=nsds/(meqn+1); 

 disp(' ')
 disp(['Frame ',num2str(Frame),' at time t = ',num2str(t)]);


 
 qmin = 1e6;
 qmax = -1e6;
 ncells = [];

 clf
 if exist('beforeframe')==2
    beforeframe  % make an m-file with this name for any other commands you
                % want executed before drawing each frame, for example
                 % if you want to use axes to specify exactly where the
                 % plot will be in the window, aspect ratio, etc.
    end 


 %=============================================
 % MAIN LOOP ON GRIDS FOR THIS FRAME:
 %=============================================
 
 % Reopen the file
 SD_id=hdfsd('start',fname,'read');
 % Initialize data set counter 
 sds_index=0; 

 for ng = 1:ngrids

   % read parameters for this grid:      
   sds_id=hdfsd('select',SD_id,sds_index);
   [data,status]=hdfsd('readdata',sds_id,0,1,17);      
   status = hdfsd('endaccess',sds_id);   
   gridno=data(1); nDim=data(2); level=data(5)+1; 
   mx = data(6); my = data(7); mz=data(8);

   xlow = data(10); ylow = data(11); zlow = data(12);
   xhigh = data(14); yhigh = data(15); zhigh = data(16);
   dx = (xhigh-xlow)/mx; dy = (yhigh-ylow)/my; dz = (zhigh-zlow)/mz;

   data=zeros([mx*my*mz,meqn]);
   start=zeros([nDim 1]); stride=ones([nDim 1]); 
   edge=stride; edge(1)=mx; edge(2)=my;
   for nq=1:meqn     
     sds_index = sds_index+1;
     sds_id=hdfsd('select',SD_id,sds_index);  
     [qcomp,status]=hdfsd('readdata',sds_id,start,stride,edge);
     status = hdfsd('endaccess',sds_id);
     data(:,nq)=reshape(qcomp,[mx*my*mz 1]);
   end      
   % Get ready for next grid   
   sds_index = sds_index+1;
   

   % if we're not plotting data at this level, skip to next grid
   if PlotData(level)==1    

   if UserVariable==1
       % User has supplied a function to convert original q variables to
       % the variable which is to be plotted, e.g. Mach number, entropy.
       qdata = feval(UserVariableFile,data);
       q = reshape(qdata,mx,my,mz);
     else
       q = reshape(data(:,mq),mx,my,mz);  
     end

   qmin = min(qmin,min(min(min(q))));  % minimum over all grids at this time
   qmax = max(qmax,max(max(max(q))));

   % keep count of how many cells at this refinement level:
   if length(ncells) < level
      ncells(level) = 0;
      end
   ncells(level) = ncells(level) + mx*my*mz;

   % -----------------------------------------------
   % plot commands go here....  First attempt with isosurface:

   hsurf = patch(isosurface(q,0.5));
   axis([0 mx 0 my 0 mz])
   set(hsurf, 'FaceColor', 'red', 'EdgeColor', 'none');
   isonormals(q,hsurf)
   camlight
   lighting phong
   grid on

   % -----------------------------------------------

   %query     % Uncomment this line to pause after plotting each subgrid
	      % Useful if you want to examine data on some subgrid

   end %   if PlotData(level)==1 
   %=============================================
   end % loop on ng (plot commands for each grid)
   %=============================================

 % done with this frame
 status = fclose(fid);

 if exist('afterframe')==2  
    afterframe  % make an m-file with this name for any other commands you
	        % want executed at the end of drawing each frame
                % for example to change the axes, or add a curve for a
                % boundary
    end

 hold off


 status = hdfsd('end',SD_id);
 end % if exist(fname)
