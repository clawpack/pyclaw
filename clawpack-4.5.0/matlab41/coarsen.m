function [qc] = coarsen(qf,ratio)
%
% function [qc] = coarsen(qf,ratio)
%
% given qf(1:mf,1:nf) on fine grid,  coarsen by a factor of ratio
% works for 1d and 2d arrays
%

[mf,nf] = size(qf);
if mf>1
    mc = mf/ratio;
    ratm = ratio;
  else
    mc = 1;
    ratm = 1;
  end
if nf>1
    nc = nf/ratio;
    ratn = ratio;
  else
    nc = 1;
    ratn = 1;
  end
indi = (0:(mc-1))*ratm + 1;
indj = (0:(nc-1))*ratn + 1;
qc = zeros(mc,nc);
for i1=0:(ratm-1)
  for j1=0:(ratn-1)
    qc = qc + qf(indi+i1, indj+j1);
    end
  end
qc = qc / (ratm*ratn);

