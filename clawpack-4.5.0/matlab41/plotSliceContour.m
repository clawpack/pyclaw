function h = plotSliceContour(xc,yc,zc,sval,sdir,xcm3,ycm3,zcm3,qcm3,clevels)

% This creates contour lines on the patch specified by (xc,yc,sval) data. We
% assume that the user has passed in (yc,zc,xval) for an xslice,
% (xc,zc,yval) for a yslice and (xc,yc,zval) for a zslice.  The contour
% levels should be a non-zero vector of values appropriate for Matlab's
% contour plotter.


% Note that we modify sval for purposes of interpolation only;  The plotting
% still happens in the true sval location.

  if (strcmp(sdir,'x'))
    % Get x slice.
    [ycm,zcm] = meshgrid(yc,zc);
    sval_new = min([max([sval,min(xc)]),max(xc)]);
    xcm = 0*ycm + sval_new;
  elseif (strcmp(sdir,'y'))
    [xcm,zcm] = meshgrid(xc,zc);
    sval_new = min([max([sval,min(yc)]),max(yc)]);
    ycm = 0*xcm + sval_new;
  elseif (strcmp(sdir,'z'))
    [xcm,ycm] = meshgrid(xc,yc);
    sval_new = min([max([sval,min(zc)]),max(zc)]);
    zcm = 0*xcm + sval_new;
  end;

  qcm = interp3(xcm3,ycm3,zcm3,qcm3,xcm,ycm,zcm,'*linear');
  if (strcmp(sdir,'x'))
    c = contourc(yc,zc,qcm,clevels,'k');
  elseif (strcmp(sdir,'y'))
    c = contourc(xc,zc,qcm,clevels,'k');
  elseif (strcmp(sdir,'z'))
    c = contourc(xc,yc,qcm,clevels,'k');
  end;

  h = drawContours(c,sval,sdir);
