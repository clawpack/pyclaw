c
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical    vtime

c
c  ::::::::::::::; ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate all grids at the input  'level' by one step of its delta(t)
c  this includes:  setting the ghost cells 
c                  advancing the solution on the grid
c                  adjusting fluxes for flux conservation step later
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      mptr = lstart(level)
      hx   = hxposs(level)
      hy   = hyposs(level)
      delt = possk(level)

 3    continue
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          locnew = node(store1,mptr)
          locaux = node(storeaux,mptr)
          time   = rnode(timemult,mptr)
c
          call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1               alloc(locaux),naux)

        mptr = node(levelptr, mptr)
        if (mptr .ne. 0) go to 3
c
c save coarse level values if there is a finer level for wave fixup
      if (level+1 .le. mxnest) then
         if (lstart(level+1) .ne. null) then
            call saveqc(level+1,nvar,naux)
         endif
      endif
c
      dtlevnew = rinfinity
      cfl_level = 0.d0    !# to keep track of max cfl seen on each level
      mptr  = lstart(level)
 5    continue
          locold = node(store2, mptr)
          locnew = node(store1, mptr)
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          time   = rnode(timemult,mptr)
c
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
c         ::: get scratch storage for fluxes and slopes
          locfp = igetsp(mitot*mjtot*nvar)
          locfm = igetsp(mitot*mjtot*nvar)
          locgp = igetsp(mitot*mjtot*nvar)
          locgm = igetsp(mitot*mjtot*nvar)
c
c  copy old soln. values into  next time step's soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
          if (level .lt. mxnest) then
             ntot   = mitot * mjtot * nvar
cdir$ ivdep
             do 10 i = 1, ntot
 10            alloc(locold + i - 1) = alloc(locnew + i - 1)
          endif
c
      xlow = rnode(cornxlo,mptr) - nghost*hx
      ylow = rnode(cornylo,mptr) - nghost*hy
      rvol = rvol + nx * ny
      rvoll(level) = rvoll(level) + nx * ny
      locaux = node(storeaux,mptr)
c
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
         locsvf = node(ffluxptr,mptr)
         locsvq = locsvf + nvar*lenbc
         locx1d = locsvq + nvar*lenbc
         call qad(alloc(locnew),mitot,mjtot,nvar,
     1            alloc(locsvf),alloc(locsvq),lenbc,
     2            intratx(level-1),intraty(level-1),hx,hy,
     3            naux,alloc(locaux),alloc(locx1d),delt,mptr)
      endif

c        # see if the grid about to advanced has gauge data to output
c        # this corresponds to previous time step, but output done
c        # now to make linear interpolation easier, since grid
c        # now has boundary conditions filled in.
c        # no testing here for mgauges>0 so that do not
c        # need to use gauges.i. the only time advanc is
c        # called that isn't "real" is in the initial setting
c        # up of grids (setgrd), but source grids are 0 there so
c        # nothing will be output.
           call dumpgauge(alloc(locnew),alloc(locaux),xlow,ylow,
     .                    nvar,mitot,mjtot,mptr)

c
      call stepgrid(alloc(locnew),alloc(locfm),alloc(locfp),
     1            alloc(locgm),alloc(locgp),
     2            mitot,mjtot,nghost,
     3            delt,dtnew,hx,hy,nvar,
     4            xlow,ylow,time,mptr,naux,alloc(locaux))

      if (node(cfluxptr,mptr) .ne. 0)
     1   call fluxsv(mptr,alloc(locfm),alloc(locfp),
     2               alloc(locgm),alloc(locgp),
     3               alloc(node(cfluxptr,mptr)),mitot,mjtot,
     4               nvar,listsp(level),delt,hx,hy)
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
         locsvf = node(ffluxptr,mptr)
         call fluxad(alloc(locfm),alloc(locfp),
     1               alloc(locgm),alloc(locgp),
     2               alloc(locsvf),mptr,mitot,mjtot,nvar,
     4               lenbc,intratx(level-1),intraty(level-1),
     5               nghost,delt,hx,hy)
      endif
c
          call reclam(locfp, mitot*mjtot*nvar)
          call reclam(locfm, mitot*mjtot*nvar)
          call reclam(locgp, mitot*mjtot*nvar)
          call reclam(locgm, mitot*mjtot*nvar)
c
          dtlevnew = dmin1(dtlevnew,dtnew)
c
          rnode(timemult,mptr)  = rnode(timemult,mptr)+delt
          mptr            = node(levelptr, mptr)
          if (mptr .ne. 0) go to 5
c
      return
      end
