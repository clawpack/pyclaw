%  infoframe1.m
%  print out info about current data

disp(['  Frame            = ' num2str(Frame)]) 
disp(['  time             = ' num2str(t)]) 
disp(['  ngrids           = ' num2str(ngrids)]) 

for ngr=1:ngrids
   disp(sprintf('  Grid %2d    level = %5d  mx = %5d  xlower = %5g  dx = %5g',  ...
             ngr,GridInfo(ngr,2), GridInfo(ngr,3),GridInfo(ngr,4), ...
	     GridInfo(ngr,5)))
   end
