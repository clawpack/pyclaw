c
c ---------------------------------------------------------
c
      logical function nestck(mnew,lbase,badpts,npts,numptc,icl,
     1                        nclust,domflags,isize,jsize,
     &      nvar,naux)
c
      implicit double precision (a-h,o-z)
      dimension  badpts(2,npts)
      integer*1  domflags(0:isize+1,0:jsize+1)

      include  "call.i"

      integer   numptc(maxcl)
c
c ::::::::::::::::::::::: NESTCK :::::::::::::::::::::::::::::::::::
c
c nestck - check that the potential grid mnew is completely
c          contained in the (coarser) finest grid which stays
c          fixed, at level lbase. projec algo. will guarantee
c          containment in all finer grids twixt them.
c          if grid not contained in some coarse grid,  then
c          bisect in long direction.
c          EVENTUALLY this has to work.
c
c input parameter:
c   mnew    - grid descriptor of potential grid
c   lbase   - level which stays fixed during regridding
c   badpts  - only the flagged pts. in this cluster (# icl)
c   domflags - the lbase domain (expanded up to level lcheck
c              which new grid must fit into
c :::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::::::::
c
       nestck = .true.
c
c      # for convex coarsest grid at level 1, nothing to check
       if (lbase .eq. 1)  go to 99

       levnew = node(nestlevel,mnew)
       lratiox  = intratx(levnew-1)
       lratioy  = intraty(levnew-1)
c
c  ### use grid indices coarsened by 1 level in checking
c  ### remember to offset by 1 since 1st grid cell is 0,0
c
       do 10 i = node(ndilo,mnew)/lratiox+1, node(ndihi,mnew)/lratiox+1
       do 10 j = node(ndjlo,mnew)/lratioy+1, node(ndjhi,mnew)/lratioy+1
          if (domflags(i,j) .eq. 0) go to 50
 10    continue
c
c      if made it here, then mnew is properly nested
       go to 99
c
c  ### grid not properly nested. bisect in long direction, and return
c  ### two clusters instead of 1.
c
 50    if (npts .gt. 1) go to 55
           write(outunit,101) levnew
           write(*,101)       levnew
 101       format(' nestck: 1 pt. cluster at level ',i5,' still not',
     1       ' nested',/,'          pt. too close to boundary')
           write(outunit,104) badpts(1,npts),badpts(2,npts)
           write(*,104)       badpts(1,npts),badpts(2,npts)
 104       format(' non-nested flagged pt. at: ',2e15.7)
           call outtre(mstart, .false.,nvar,naux)
           call outmsh(mnew, .false.,nvar,naux)
           stop

 55    if (nclust .lt. maxcl) go to 60
           write(outunit,102) maxcl
           write(*,102)       maxcl
 102       format(' too many clusters: > ',i5,' (from nestck) ')
           stop

 60   if (nprint) write(outunit,103) icl, npts
 103  format(' bisecting cluster ',i5,' with ',i5,' pts. in nestck')
      if (rnode(cornxhi,mnew)-rnode(cornxlo,mnew) .gt.
     1    rnode(cornyhi,mnew) - rnode(cornylo,mnew)) then
           rmid  = (rnode(cornxhi,mnew) + rnode(cornxlo,mnew) ) / 2.
           rmid  = (node(ndihi,mnew) + node(ndilo,mnew) + 1 ) / 2.
           rmid  = rmid / lratiox
           idir = 1
       else
           rmid  = (rnode(cornyhi,mnew) + rnode(cornylo,mnew) ) / 2.
           rmid  = (node(ndjhi,mnew) + node(ndjlo,mnew) + 1) / 2.
           rmid  = rmid / lratioy
           idir = 2
       endif
c
       ipt  = 1
       ntop = npts

 90    if (badpts(idir,ipt) .lt. rmid) go to 100
c
c  ### swap with a point in top half not yet tested. keep smaller
c  ### half of rect. in bottom half
c
       temp           = badpts(1,ipt)
       badpts(1,ipt)  = badpts(1,ntop)
       badpts(1,ntop) = temp
       temp           = badpts(2,ipt)
       badpts(2,ipt)  = badpts(2,ntop)
       badpts(2,ntop) = temp
       ntop           = ntop - 1
       if (ipt .le. ntop) go to 90
       go to 110
 100   ipt = ipt +1
       if (ipt .le. ntop) go to 90
c
c  ### ntop points to top of 1st cluster (= no. of points in 1st cluster)
c
 110   numptc(icl)     = npts - ntop
       do 120 i        = icl, nclust
          nmove           = nclust + icl - i
 120      numptc(nmove+1) = numptc(nmove)
       numptc(icl)     = ntop
       nclust          = nclust + 1
       nestck          = .false.
c
 99    return
       end
