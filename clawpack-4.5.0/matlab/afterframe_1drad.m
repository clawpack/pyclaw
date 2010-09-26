% Sample afterframe.m file for 2d or 3d that compares scatter plot results
% plotted using PlotType=4 to a fine-grid 1d result in the subdirectory
% 1dref (e.g. for radially-symmetric problems).   A function map1d.m
% may also need to be provided so that the scatter plot uses the appropriate
% 1d variable.
%

if PlotType==4
   dir = './1dref/';
   dim = 1;
   [amrdata1d,t1d] = readamrdata(dim,Frame,dir);
   if isempty(t1d)
      disp('Run xclaw in 1dref to generate 1d reference solution')
    else
      hold on;
      [q1d,x1d] = plotframe1ez(amrdata1d,mq,'r-');
      hold off;
    end
  end

