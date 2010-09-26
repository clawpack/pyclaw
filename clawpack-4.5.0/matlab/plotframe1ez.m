function [q,xc,h]  = plotframe1ez(amrdata,mq,plotstyle,...
    uservariablefile,mapc2pfile)

% PLOTFRAME1EZ for plotting simple 1d plots
%
%   PLOTFRAME1EZ is a functional version of plotframe1 and can be used in
%   conjunction with other 2d and 3d plotting routines without overwriting
%   global variables from those routines.
%
%   [Q,X] = PLOTFRAME1EZ(AMRDATA,MQ,PLOTSTYLE,USERVARIABLEFILE,MAPC2PFILE)
%   plots the data stored in AMRDATA, using the variables MQ,PLOTSTYLE,
%   USERVARIABLEFILE and MAPC2PFILE to specify the type of plot to show.
%   AMRDATA is created with the READAMRDATA command, and is a structure
%   containing grid information.  See READAMRDATA.  MQ is a scalar which
%   specifies which variable (i.e. column) of the amrdata.data matrix to
%   read.  USERVARIABLEFILE is a string specifing any user defined function
%   which the data is passed through.  It has the syntax :
%
%             data = amrdata.data';
%             q = feval(UserVariableFile,data);
%
%   If USERVARIABLEFILE is the empty string ('') then q is set to data(:,mq).
%
%   Finally, if the user has set mappedgrid = 1, then a function 'mapc2p' is
%   called to map the Cartesian data to physical data.  MAPC2P should
%   be the name of a function on the Matlab path with the syntax
%
%         xp = mapc2p(xc)
%
%    NOTE : It is not assumed that this function has the name 'mapc2p',
%    because this would conflict with the mapc2p used for the 2d or 3d
%    routines.
%
%    PLOTFRAME1EZ returns the 1d Q and X data.
%
%   [Q,X] = PLOTFRAME1EZ(AMRDATA,MQ,PLOTSTYLE) uses default values
%   ('') for USERVARIABLEFILE and MAPC2PFILE.
%
%   [Q,X] = PLOTFRAME1EZ(AMRDATA,MQ) sets PLOTSTYLE = 'b';
%
%   [Q,X] = PLOTFRAME1EZ(AMRDATA) sets MQ = 1.
%
%   [Q,X,P] = PLOTFRAME1EZ(....) returns a handle P to the 1d line plot.
%
%    Example :
%
%        % To call this routine from an afterframe.m file
%        % ... <afterframe>...
%        hold on;
%        dir = './results1d/';
%        dim = 1;
%        [amrdata1d,t1d] = readamrdata(dim,Frame,dir);
%        userfile = '';  % UserVariableFile
%        mapc2pfile = 'mapc2p1d';
%        [q1d,x1d] = plotframe1ez(amrdata1d,mq,'b-',userfile,mapc2pfile);
%        hold off;
%
%    See also SETPLOT, MappedGrid, READAMRDATA.

if (length(amrdata) == 0)
  disp('*** plotframe1ez : ''amrdata'' is empty');
  return;
end;

if (length(amrdata) > 1)
  disp('*** plotframe1ez : length(amrdata) > 1; plotframe1ez does not ');
  disp('*** support multiple grids');
  return;
end;

if (nargin < 5)
  mapc2pfile = '';
  if (nargin < 4)
    uservariablefile = '';
    if (nargin < 3)
      plotstyle = 'r-';
      if (nargin < 2)
	mq = 1;
      end;
    end;
  end;
end;

uservariable = ~isempty(uservariablefile);
mappedgrid = ~isempty(mapc2pfile);

gridno = amrdata(1).gridno;
level = amrdata(1).level;
mx = amrdata(1).mx;
xlow = amrdata(1).xlow;
dx = amrdata(1).dx;
data = amrdata(1).data';

xc = xlow + ((1:mx) - 0.5)*dx;

% for compatibility with old matlab41/plotframe1 convention:
x = xc;

if (uservariable == 1)
  % User has supplied a function to convert original q variables to
  % the variable which is to be plotted, e.g. Mach number, entropy.
  q = feval(uservariablefile,data);
else
  q = data(:,mq);
end;

if mappedgrid
  % coordinate mapping must be applied
  xc = feval(mapc2pfile,xc);
end

[lstyle,lcolors,mstyle] = get_plotstyle({plotstyle},1);

p = line('XData',xc,'YData',q,'LineStyle',lstyle{1},...
    'Color',lcolors{1},'Marker',mstyle{1});

if (nargout == 3)
  % sent out plot handle.
  h = p;
end;
