%  infoplot2.m
%  print out parameter values

disp(['  PlotType         = ' num2str(PlotType)])
disp(['  mq               = ' num2str(mq)])
disp(['  UserVariable     = ' num2str(UserVariable)])
disp(['  UserVariableFile = ' UserVariableFile])
disp(['  MappedGrid       = ' num2str(MappedGrid )])
disp(['  PlotData         = ' num2str(PlotData)])
disp(['  PlotGrid         = ' num2str(PlotGrid)])
disp(['  PlotGridEdges    = ' num2str(PlotGridEdges)])


if (exist('ContourValues'))
  disp(['  ContourValues    = ' num2str(ContourValues)]);
end;

if (exist('IsosurfValues'))
  disp(['  IsosurfValues    = ' num2str(IsosurfValues)]);
end;
