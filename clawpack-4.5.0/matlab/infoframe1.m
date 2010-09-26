%  infoframe1.m
%  print out info about current data

disp(['  Frame            = ' num2str(Frame)])
disp(['  time             = ' num2str(t)])
disp(['  ngrids           = ' num2str(ngrids)])

for ngr=1:ngrids
   disp(sprintf('  Grid %2d    level = %5d  mx = %5d  xlower = %5g  dx = %5g',  ...
             ngr,amrdata(ngr).level, amrdata(ngr).mx, amrdata(ngr).xlower, ...
	     amrdata(ngr).dx);
   end
