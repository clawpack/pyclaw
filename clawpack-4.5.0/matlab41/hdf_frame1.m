% claw/matlab/plotframe1.m
%
% Plot data at a single time produced by clawpack, 
%     files at this time frame
% Assumes Frame is already set.  
%
% Assumes variables described in claw/matlab/setplot1.m are set, e.g. 
% by calling that routine.

% read time and number of grids:
n1 = Frame+10000;
fname = ['q',num2str(n1),'.hdf'];
fname(2) = '0';
 
if ~exist(fname) 
     disp(' ');
     disp(['Frame ',num2str(Frame),' does not exist ***']);
     end

if exist(fname)

 % start reading data, beginning with parameters for this frame:

 SD_id=hdfsd('start',fname,'read');
 [nsds,nattr,status]=hdfsd('fileinfo',SD_id); 
 sds_index=0;
 sds_id=hdfsd('select',SD_id,sds_index);
 [data,status]=hdfsd('readdata',sds_id,0,1,16);
 status = hdfsd('endaccess',sds_id);
 status = hdfsd('end',SD_id);
 t = data(3);
 meqn = data(4);
 ngrids=nsds/(meqn+1); 

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
   [data,status]=hdfsd('readdata',sds_id,0,1,16);
   status = hdfsd('endaccess',sds_id);     
   gridno=data(1); nDim=data(2); level=data(5)+1; 
   mx = data(6); 

   xlow = data(10);
   xhigh = data(14);
   dx = (xhigh-xlow)/mx;

   GridInfo(ng,1) = gridno;
   GridInfo(ng,2) = level;
   GridInfo(ng,3) = mx;
   GridInfo(ng,4) = xlow;
   GridInfo(ng,5) = dx;

   % read q data:   
   
   data=zeros([mx,meqn]);
   start=zeros([nDim 1]); stride=ones([nDim 1]); 
   edge=stride; edge(1)=mx; 
   for nq=1:meqn
     sds_index = sds_index+1;     
     sds_id=hdfsd('select',SD_id,sds_index);  
     status = hdfsd('endaccess',sds_id);
     [qcomp,status]=hdfsd('readdata',sds_id,start,stride,edge);
     data(:,nq)=reshape(qcomp,[mx 1]);
   end   
   
   % Get ready for next grid   
   sds_index = sds_index+1;


   if UserVariable==1
       % User has supplied a function to convert original q variables to
       % the variable which is to be plotted, e.g. Mach number, entropy.
       q = feval(UserVariableFile,data);
       end

   x = xlow + ((1:mx) - 0.5)*dx;

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
  else if length(mq)==1
     % plotting a single variable
     q = data(:,mq);
     plot(x,q,plotstyle)
     title(['q(',num2str(mq),')   at time t = ', num2str(t)])
  else
     % plotting more than one variable in same figure
     % use subplot to put each in
     for  mq1=mq
        subplot(length(mq),1,find(mq==mq1))
        q = data(:,mq1);
        plot(x,q,plotstyle)
        title(['q(',num2str(mq1),')   at time t = ', num2str(t)])
        end  % loop on mq1
     end  % if UserVariable
  end  % loop on ng
 

 if exist('afterframe')==2  
    hold on
    afterframe  % make an m-file with this name for any other commands you
	        % want executed at the end of drawing each frame
                % for example to change the axes, or add something to the plot
    hold off
    end


 end
 
 status = hdfsd('end',SD_id);
end % if exist(fname)
