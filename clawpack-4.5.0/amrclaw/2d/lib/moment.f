c
c ----------------------------------------------------------
c
      subroutine moment (intrect,badpts,npt,usage)
c
      implicit double precision (a-h,o-z)

      include "call.i"

      dimension     intrect(nsize),badpts(2,npt)
c
c :::::::::::::::::::::::: MOMENT ::::::::::::::::::::::::::::::::::
c  moment = compute enclosing rectangle around flagged points.
c  save some info., even tho. usage might be low and rect. scrapped.

c input parameters:
c     badpts      = x,y coords of flagged badpts grouped into clusters
c                   are in the first two rows
c     npt         = num. of badpts. in the cluster.
c
c output parameters:
c     usage       = ratio of flagged to unflagged badpts. in new grid
c                   measures goodness of fit and clustering
c    intrect( )    = stores some info. for grid created herein.
c                   sometimes rect = rnode, sometimes = temp. array.
c                   sometimes intrect = node.
c                   depending on calling prog. (grdfit or expand)
c
c
c :::::::::::::::::::::::: MOMENT ::::::::::::::::::::::::::::::::::
c
      rn = dble(npt)
c
c compute length of enclosing rectangles to include all flagged badpts.
c
      emx1 = badpts(1,1)
      emn1 = emx1
      emx2 = badpts(2,1)
      emn2 = emx2
      do 80 ipt = 1, npt
          if (badpts(1,ipt) .gt. emx1) emx1 = badpts(1,ipt)
          if (badpts(1,ipt) .lt. emn1) emn1 = badpts(1,ipt)
          if (badpts(2,ipt) .gt. emx2) emx2 = badpts(2,ipt)
          if (badpts(2,ipt) .lt. emn2) emn2 = badpts(2,ipt)
 80   continue
c
c from length of the sides determine rect. corners.
c transform to cell numbers (subtract .5)
c
      intrect(ndilo) = nint(emn1 - .5)
      intrect(ndjlo) = nint(emn2 - .5)
      intrect(ndihi) = nint(emx1 - .5)
      intrect(ndjhi) = nint(emx2 - .5)
c
c compute usage
c
      iside1 = intrect(ndihi) - intrect(ndilo) + 1
      iside2 = intrect(ndjhi) - intrect(ndjlo) + 1
      gpall  = iside1 * iside2
      usage  = rn / gpall
c
      return
      end
