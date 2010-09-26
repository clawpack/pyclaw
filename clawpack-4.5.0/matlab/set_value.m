function set_value(a,a_user,a_default)

% Internal Matlab routine for Clawpack graphics

% This routine executes the following code in the base
% workspace :
%
%   if (~exist(a_user))
%      a = a_default;
%   else
%      eval(['a = ',a_user]);
%   end;
%
%    Here, a and a_user are the names of variables;  a_default is a value.
%
% Example :
%
%           set_value(mappedgrid,'MappedGrid',0);
%
% This is equivilent to
%
%      if (exist('MappedGrid'))
%          mappedgrid = MappedGrid;
%      else
%          mappedgrid = 0;
%      end;
%

if (~evalin('base',['exist(''',a_user,''')']))
  assignin('base',a,a_default);
else
  evalin('base',[a, ' = ',a_user,';']);
end;
