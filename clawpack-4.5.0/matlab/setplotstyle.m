function syms = setplotstyle(varargin)

% SETPLOTSTYLE sets the graphics symbols to be used for 1d plots and
%    2d/3d scatter plots.
%
%    S = SETPLOTSTYLE(s1,s2,s3...) takes as input arguments potential strings
%    to the PLOT command and returns a cell matrix S to be used by the
%    functions PLOTCLAW1, or PLOTCLAW2 and PLOTCLAW3 when PlotType = 4 is set.
%    Number of input arguments corresponds to number of desired AMR levels
%    to plot.
%
%    Example :
%
%        % For 1d plots, set PlotStyle for three levels of amr data.
%        PlotStyle = setplotstyle('r-','y*','g--');
%
%    Example :
%
%        % For 2d/3d scatter plots or line plots, set ScatterStyle
%        % (or, equivalently, LineStyle).for three levels of amr data
%        ScatterStyle = setplotstyle('b*','gx','r.');
%
%     This assures that the variables is set properly for use
%     with the graphics routines, and is equivalent setting ScatterStyle
%     (or PlotStyle or LineStyle) to a cell matrix :
%
%         ScatterStyle = {'b*','r--','gx','r^'};
%
%     See also PLOTCLAW1, PLOTCLAW2, PLOTCLAW3, SETPLOT.


n = nargin;
syms = cell(n,1);

for i = 1:n,
  syms{i} = varargin{i};
end;
