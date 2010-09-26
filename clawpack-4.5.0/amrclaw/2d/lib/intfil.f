c
c ------------------------------------------------------
c
       subroutine intfil(val,mi,mj,time,flaguse,nrowst,ncolst,
     2                   ilo,ihi,jlo,jhi,level,nvar,naux)
c
c ::::::::::::::::::::::: INTFIL ::::::::::::::::::::::::::::::::;
c  INTFIL: interpolates values for a patch at the specified level and
c  location, using values from grids at LEVEL and coarser, if nec.
c
c  take the intersection of a grid patch with corners at ilo,ihi,jlo,jhi
c  and all grids mptr at LEVEL.  If there is a non-null intersection
c  copy the solution vaues from mptr (at TIME) into VAL array.
c  assumes patch at same level so do straight copy, not skipping
c  every intrat or doing any interpolation here,
c  assume called in correct order of levels, so that when copying
c  is ok to overwrite.
c
c  N.B.: there are no dummy points around patch, since
c        this is not an official "grid" but a "patch".
c
c  used array marks when point filled. filpatch checks if points left over
c  after intersections at specified level.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
      dimension   val(mi,mj,nvar)
      logical     tinterp
c
      iadd(i,j,ivar)   = loc    + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iadnew(i,j,ivar) = locnew + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iadold(i,j,ivar) = locold + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      dimension flaguse(ilo:ihi, jlo:jhi)
c
      dt     = possk(level)
c     teps   = dt / 10.d0
c     need non-dimensional epsilon for time. was a problem
c     in large scaled tsunami tests
      teps   = 1.d-4



      do i = ilo, ihi
      do j = jlo, jhi
        flaguse(i,j) = 0.0
      end do
      end do

      mptr   = lstart(level)
 10   if (mptr .eq. 0) go to 105
c
c     :::  check if grid mptr and patch intersect
c
      imlo = node(ndilo, mptr)
      jmlo = node(ndjlo, mptr)
      imhi = node(ndihi, mptr)
      jmhi = node(ndjhi, mptr)

      nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      mitot = nx + 2*nghost
      mjtot = ny + 2*nghost

      ixlo = max(imlo,ilo)
      ixhi = min(imhi,ihi)
      jxlo = max(jmlo,jlo)
      jxhi = min(jmhi,jhi)
c
      if (.not.(ixlo .le. ixhi .and. jxlo .le. jxhi)) go to 90
c
c  :::  grids intersect. figure out what time to use.
c  :::  alphai = 1 for new time; 0 for old time
c
      alphac = (rnode(timemult,mptr) - time)/dt
      alphai = 1.d0-alphac

      if ((alphai .lt. -teps) .or. (alphai .gt. 1.d0+teps)) then
          write(outunit,900) time, mptr, level
          write(*,900) time, mptr, level
 900      format(' time wanted ',e15.7,' not available from grid ',i4,
     1           'level',i4)
          write(outunit,'(A,E24.16)') 'Line 80', dt
          write(outunit,901) ilo,ihi,jlo,jhi,mptr,level,time,
     .                 rnode(timemult,mptr),alphai,teps
          write(*,901) ilo,ihi,jlo,jhi,mptr,level,time,
     .                 rnode(timemult,mptr),alphai,teps
          call outtre(mstart,.false.,nvar,naux)
          stop
      endif
c
      tinterp = .false.
      if (dabs(alphai - 1.d0) .lt. teps) then
          loc = node(store1,mptr)
      else if (dabs(alphai) .lt. teps) then
          loc = node(store2,mptr)
          if (level .eq. mxnest) then
             write(outunit,'(A,E24.16)') 'Line 95', dt
             write(outunit,901) ilo,ihi,jlo,jhi,mptr,level,time,
     .                          rnode(timemult,mptr),alphai,teps
             write(*,901) ilo,ihi,jlo,jhi,mptr,level,time,
     .                          rnode(timemult,mptr),alphai,teps
             stop
           endif
      else
          locold  = node(store2,mptr)
          locnew  = node(store1,mptr)
          tinterp = .true.
          if (level .eq. mxnest) then
             write(outunit,'(A,E24.16)') 'Line 107',dt
             write(outunit,901) ilo,ihi,jlo,jhi,mptr,level,time,
     .                    rnode(timemult,mptr),alphai,teps
             write(*,901) ilo,ihi,jlo,jhi,mptr,level,time,
     .                    rnode(timemult,mptr),alphai,teps
             stop
          endif
      endif
 901  format(' trying to interpolate from previous time values ',/,
     .       ' for a patch with corners ilo,ihi,jlo,jhi:'
     .       ,/,2x,4i10,/,
     .       ' from source grid ',i4,' at level ',i4,/,
     .       ' time wanted ',e24.16,' source time is ',e24.16,/,
     .       ' alphai, teps ',2e24.16)
c
      if (.not. tinterp) then
c     ::: no time interp. copy the solution values
         do 45 ivar = 1, nvar
         do 35 j = jxlo, jxhi
         do 20 i = ixlo, ixhi
             val(i-ilo+nrowst,j-jlo+ncolst,ivar) =
     1            alloc(iadd(i-imlo+nghost+1,j-jmlo+nghost+1, ivar))
             flaguse(i,j) = 1.d0
 20      continue
 35      continue
 45      continue
      else
c     ::: linear interpolation in time
         do 85 ivar = 1, nvar
         do 75 j = jxlo, jxhi
         do 65 i = ixlo, ixhi
           val(i-ilo+nrowst,j-jlo+ncolst,ivar) =
     1      alloc(iadnew(i-imlo+nghost+1,j-jmlo+nghost+1,ivar))*alphai +
     2      alloc(iadold(i-imlo+nghost+1,j-jmlo+nghost+1,ivar))*alphac
            flaguse(i,j) = 1.d0
 65      continue
 75      continue
 85      continue
      endif
c
 90   mptr = node(levelptr, mptr)
      go to 10
c
 105  continue

c  set used array points which intersect boundary to be equal to 1;
c  they will be set elsewhere

      if (jhi .ge. jregsz(level)) then
        do 1000 j = max(jregsz(level),jlo), jhi
        do 1000 i = ilo, ihi
           flaguse(i,j) = 1.d0
1000    continue
      endif

      if (jlo .lt. 0) then
        ncolend = ncolst + jhi - jlo 
        do 1200 j = jlo, min(-1,ncolend)
        do 1200 i = ilo, ihi
           flaguse(i,j) = 1.d0
1200    continue
      endif

      if (ilo .lt. 0) then
        nrowend = nrowst + ihi - ilo
        do 1400 i = ilo, min(-1,nrowend)
        do 1400 j = jlo, jhi
           flaguse(i,j) = 1.d0
1400    continue
      endif

      if (ihi .ge. iregsz(level)) then
        do 1600 i = max(iregsz(level),ilo), ihi
        do 1600 j = jlo, jhi
           flaguse(i,j) = 1.d0
1600    continue
      endif
c
      return
      end
