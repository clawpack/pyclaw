%  infoplot1.m
%  print out parameter values

disp(['  mq               = ' num2str(mq)])
disp(['  UserVariable     = ' num2str(UserVariable)])
disp(['  UserVariableFile = ' UserVariableFile])
disp(['  MappedGrid       = ' num2str(MappedGrid )])

if (exist('plotstyle'))
  disp(['  plotstyle        = ' plotstyle])
end;
if (exist('PlotStyle'))
  disp(['  PlotStyle        = ' PlotStyle]);
end;
