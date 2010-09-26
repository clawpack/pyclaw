function bn = getblocknumber()

% GETBLOCKNUMBER returns current block number.
%
%      BN = GETBLOCKNUMBER returns the block number (or patch) currently being
%      created.  This is used primarly for multi-block calculations in which
%      a mapped grid is being used.  Typically, the mapping used to map from
%      Cartesian coordinates to physical coordinates (the 'mapc2p' function)
%      will depend on the block number.
%
%      This will usually be called within a 'mapc2p' routine.
%
% See also mappedgrid.
%

% The block number is set by set_blocknumber, which is called from
% plotframe2, ploframe3 and showmesh.m

global bn_set_blocknumber;

bn = bn_set_blocknumber;
