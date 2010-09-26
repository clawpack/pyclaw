function [amr,t] = readamrdata_hdf(dim,Frame,dir);

% Internal routines for Clawpack graphics.

% Check to see if an hdf file exists for this frame.
n1 = Frame+10000;
fname = [dir, 'fort.',num2str(n1),'.hdf'];
fname(length(dir)+6) = 'q';

if ~exist(fname)
  amr = {};
  t = [];
  return;
end

% Read grid info for initial grid from hdf file.
disp(['Reading data from ',fname]);
sd_id = hdfsd('start',fname,'read');

% Get number of datasets (nsds) and global attributes (nattr).
[nsds,nattr,status] = hdfsd('fileinfo',sd_id);

% Select and read in the first dataset (grid info for grid 0).
sds_index = 0;
sds_id = hdfsd('select',sd_id,sds_index);
[data,status] = hdfsd('readdata',sds_id,0,1,21);
status = hdfsd('endaccess',sds_id);

% Extract time, meqn and ngrids.
t      = data(3);
meqn   = data(4);
ngrids = nsds/(meqn+1);

%=============================================
% MAIN LOOP ON GRIDS FOR THIS FRAME:
%=============================================

% Initialize data set counter
sds_index=0;

for ng = 1:ngrids

  % get parameters for this grid from hdf file.
  sds_id = hdfsd('select',sd_id,sds_index);
  [data,status] = hdfsd('readdata',sds_id,0,1,21);
  status = hdfsd('endaccess',sds_id);

  % check number of dimensions requested matches that in HDF file.
  if dim ~= data(2)
    error(['Number of dimensions requested does not match HDF ' ...
           'file']);
  end

  % copy parameters for this grid into amrdata data structure.
  % note that data(14:17) contains xhigh, yhigh, zhigh.
  amrdata.gridno = data(1);
  amrdata.level  = data(5);

  if (dim == 1)
    amrdata.mx   = data(6);
    amrdata.xlow = data(10);
    amrdata.dx   = data(18);

    % create an array of zeros which will hold the data.
    amrdata.data = zeros([meqn,amrdata.mx]);

  elseif (dim == 2)
    amrdata.mx   = data(6);
    amrdata.xlow = data(10);
    amrdata.dx   = data(18);

    amrdata.my   = data(7);
    amrdata.ylow = data(11);
    amrdata.dy   = data(19);

    % create an array of zeros which will hold the data.
    amrdata.data = zeros([meqn,amrdata.mx*amrdata.my]);

  elseif (dim == 3)

    amrdata.mx   = data(6);
    amrdata.xlow = data(10);
    amrdata.dx   = data(18);

    amrdata.my   = data(7);
    amrdata.ylow = data(11);
    amrdata.dy   = data(19);

    amrdata.mz   = data(8);
    amrdata.zlow = data(12);
    amrdata.dz   = data(20);

    % create an array of zeros which will hold the data.
    amrdata.data=zeros([meqn,amrdata.mx*amrdata.my*amrdata.mz]);

  end

  % provide description of the data set for reading data.
  start  = zeros([dim,1]);
  stride = ones([dim,1]);
  if (dim == 1)
    edge = [amrdata.mx]';
  elseif (dim == 2)
    edge = [amrdata.mx amrdata.my]';
  elseif (dim == 3)
    edge = [amrdata.mx amrdata.my amrdata.mz]';
  end

  % loop over datasets on each grid, reading each one.
  for nq = 1:meqn
    sds_index = sds_index+1;
    sds_id         = hdfsd('select',sd_id,sds_index);
    [qcomp,status] = hdfsd('readdata',sds_id,start,stride,edge);
    status         = hdfsd('endaccess',sds_id);
    if (dim == 1)
      amrdata.data(nq,:) = reshape(qcomp,[amrdata.mx 1])';
    elseif (dim == 2)
      amrdata.data(nq,:) = reshape(permute(qcomp,[2 1]),...
                                   [amrdata.mx*amrdata.my 1])';
    elseif (dim == 3)
      amrdata.data(nq,:) = reshape(permute(qcomp,[3 2 1]),...
                                   [amrdata.mx*amrdata.my*amrdata.mz 1])';
    end

  end

  % Get ready for next grid
  sds_index = sds_index+1;

  amr(ng) = amrdata;

end

status = hdfsd('end',sd_id);
