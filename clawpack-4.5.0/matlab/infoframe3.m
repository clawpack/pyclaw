
%  infoframe2.m
%  print out info about current data

disp(['  Frame            = ' num2str(Frame)])
disp(['  time             = ' num2str(t)])

for ngr=1:ngrids
   disp(sprintf('  Grid %2d    level = %5d  mx = %5d  xlower = %5g  dx = %5g',  ...
             ngr,amrdata(ngr).level,amrdata(ngr).mx, amrdata(ngr).xlow,...
	     amrdata(ngr).dx));
   disp(sprintf('                            my = %5d  ylower = %5g  dy = %5g',  ...
             amrdata(ngr).my, amrdata(ngr).ylow,amrdata(ngr).dy));
   disp(sprintf('                            mz = %5d  zlower = %5g  dz = %5g',  ...
             amrdata(ngr).mz, amrdata(ngr).zlow,amrdata(ngr).dz));
   end


if ngrids>1
  disp(['  number of grid cells on each level:  ',num2str(ncells)]);
  end
disp(['  total number of grid cells:  ',num2str(sum(ncells))]);
disp(['  q values range between:  ',num2str(qmin),' and ',num2str(qmax)]);
disp(' ')
