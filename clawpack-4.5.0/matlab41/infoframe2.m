
%  infoframe2.m
%  print out info about current data

disp(['  Frame            = ' num2str(Frame)]) 
disp(['  time             = ' num2str(t)]) 

for ngr=1:ngrids
   disp(sprintf('  Grid %2d    level = %5d  mx = %5d  xlower = %5g  dx = %5g',  ...
             GridInfo(ngr,1),GridInfo(ngr,2), GridInfo(ngr,3),GridInfo(ngr,5), ...
	     GridInfo(ngr,7)))
   disp(sprintf('                            my = %5d  ylower = %5g  dy = %5g',  ...
             GridInfo(ngr,4), GridInfo(ngr,6),GridInfo(ngr,8)))
   end


if ngrids>1
  disp(['  number of grid cells on each level:  ',num2str(ncells)]);
  end
disp(['  total number of grid cells:  ',num2str(sum(ncells))]);
disp(['  q values range between:  ',num2str(qmin),' and ',num2str(qmax)]);
disp(' ')
