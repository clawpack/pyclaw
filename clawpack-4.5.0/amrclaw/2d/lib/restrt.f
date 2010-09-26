
c
c ---------------------------------------------------------
c
      subroutine restrt(nsteps,time,nvar)
c
      implicit double precision (a-h,o-z)
      logical   ee
 
      include  "call.i"
 
      character*12 rstfile
      logical foundFile
      dimension intrtx(maxlv),intrty(maxlv),intrtt(maxlv)
c
c :::::::::::::::::::::::::::: RESTRT ::::::::::::::::::::::::::::::::
c read back in the check point files written by subr. check.
c
c some input variables might have changed, and also the
c alloc array could have been written with a smaller size at checkpoint
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      rstfile  = 'restart.data'
      inquire(file=rstfile,exist=foundFile)
      if (.not. foundFile) then
        write(*,*)" Did you forget to link file 'restart.data'"
        write(*,*)"  to a checkpoint file?"
        stop
      endif
      open(rstunit,file=rstfile,status='old',form='unformatted')
      rewind rstunit

      read(rstunit) lenmax,lendim,isize

c     # need to allocate for dynamic memory:
      call restrt_alloc(isize)

      read(rstunit) (alloc(i),i=1,lendim)
      read(rstunit) hxposs,hyposs,possk,icheck
      read(rstunit) lfree,lenf
      read(rstunit) rnode,node,lstart,newstl,listsp,tl,
     1       ibuf,mstart,ndfree,lfine,iorder,mxnold,
     2       intrtx,intrty,intrtt,iregsz,jregsz,
     2       iregst,jregst,iregend,jregend,
     3       kcheck1,nsteps,time,
     3       matlabu
      read(rstunit) evol, rvol, rvoll, lentot

      write(outunit,100) nsteps,time
      write(6,100) nsteps,time
 100  format(/,' RESTARTING the calculation after ',i5,' steps',
     1        /,'  (time = ',e15.7,')')
c
c     error checking that refinement ratios have not changed
c
      do i = 1, mxnold-1
        if ( (intratx(i) .ne. intrtx(i)) .or.
     .       (intraty(i) .ne. intrty(i)) .or.
     .       (kratio(i) .ne.  intrtt(i)) ) then
        write(outunit,*) 
     .  " not allowed to change existing refinement ratios on Restart"
        write(*,*)
     .  " not allowed to change existing refinement ratios on Restart"
        write(outunit,*)" Old ratios:"
        write(*,*)      " Old ratios:"
        write(outunit,903)(intrtx(j),j=1,mxnold-1)
        write(*,903)      (intrtx(j),j=1,mxnold-1)
        write(outunit,903)(intrty(j),j=1,mxnold-1)
        write(*,903)      (intrty(j),j=1,mxnold-1)
        write(outunit,903)(intrtt(j),j=1,mxnold-1)
        write(*,903)      (intrtt(j),j=1,mxnold-1)
 903    format(6i3)
        stop
       endif
      end do

c
c adjust free list of storage in case size has changed.
c
      idif = memsize - isize
      if (idif .gt. 0) then
          lfree(lenf,1) = isize + 2
          call reclam(isize+1,idif)
      else if (idif .lt. 0) then
            write(outunit,900) isize, memsize
            write(*,900)       isize, memsize
 900        format(' size of alloc not allowed to shrink with ',/,
     .             ' restart old size ',i7,' current size  ',i7)
            stop
      endif
c
c adjust storage in case mxnest has changed - only allow it to increase,
c
       if (mxnest .eq. mxnold) go to 99

       if (mxnest .lt. mxnold) then
         if (lfine .lt. mxnest) then
             go to 99
         else
             write(outunit,901) mxnold, mxnest
             write(*,      901) mxnold, mxnest
901          format(' only allow mxnest to increase: ',/,
     &            '  old mxnest ',i4, ' new mxnest ',i4)
              stop
         endif
       endif

c add second storage loc to grids at previous mxnest
        ee = .false.
        do 10 level = 1, mxnold
           if (icheck(level) .ge. kcheck) then
              ee = .true.
           endif
           write(*,*)" increasing max num levels from ",mxnold,
     .                 ' to',mxnest
           write(outunit,*)" increasing max num levels from ",mxnold,
     .                 ' to',mxnest
10      continue

        if (ee .and. tol .gt. 0) then
c           ## if Richardson used, will delay error est. 1 step til have old soln. vals
             write(*,*)" first Richardson error estimation step"
             write(*,*)" will estimate mostly spatial error "
             write(outunit,*)" first Richardson error estimation step"
             write(outunit,*)" will estimate mostly spatial error  "
         endif

c          #  add second storage location to previous mxnest level
         mptr = lstart(mxnold)
15       if (mptr .eq. 0) go to 25
            mitot = node(ndihi,mptr)-node(ndilo,mptr)+1+2*nghost
            mjtot = node(ndjhi,mptr)-node(ndjlo,mptr)+1+2*nghost
            node(store2,mptr) = igetsp(mitot*mjtot*nvar)
            mptr = node(levelptr,mptr)
            go to 15
25       continue
c
c          # add new info. to spatial and counting arrays
 99        level = lfine + 1
           rrk = dble(kratio(lfine))
35         if (level .gt. mxnest) go to 45
             hxposs(level) = hxposs(level-1) / dble(intratx(level-1))
             hyposs(level) = hyposs(level-1) / dble(intraty(level-1))
             possk (level) = possk (level-1) / rrk
             iregsz(level) = iregsz(level-1) * intratx(level-1)
             jregsz(level) = jregsz(level-1) * intraty(level-1)
             rrk           = kratio(level)
             level         = level + 1
             go to 35
45         continue
c
c
      return
      end
