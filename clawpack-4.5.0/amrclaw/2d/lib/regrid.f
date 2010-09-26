c
c -----------------------------------------------------------
c
      subroutine regrid  (nvar,lbase,cut,naux,t0)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
c :::::::::::::::::::::::::::: REGRID :::::::::::::::::::::::::::::::

c  regrid = flag points on each grid with a level > = lbase.
c  cluster them, and fit new subgrids around the clusters.
c  the lbase grids stay fixed during regridding operation.
c  when a parent grid has its error estimated, add its kid grid
c  information to the error grid before clustering. (project)
c  order of grid examination - all grids at the same level, then
c  do the next coarser level.
c
c input parameters:
c     lbase  = highest level that stays fixed during regridding
c     cutoff = criteria for measuring goodness of rect. fit.

c local variables:
c     lcheck = the level being examined.
c     lfnew  = finest grid to be. will replace lfine.

c global
c    mstart  = start of very coarsest grids.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      lcheck    = min0(lfine,mxnest-1)
      lfnew     = lbase
      do 10 i   = 1, mxnest
 10     newstl(i) = 0
      time      = rnode(timemult, lstart(lbase))
c
 20   if (lcheck .lt. lbase) go to 50
          call grdfit(lbase,lcheck,nvar,naux,cut,time,t0)
          if (newstl(lcheck+1) .eq. 0) go to 40
          lfnew = max0(lcheck + 1,lfnew)
 40       continue
          lcheck = lcheck - 1
c
      go to 20
 50   continue
c
c  end of level loop
c
c  remaining tasks left in regridding:
c
c  interpolate storage for the new grids.  the starting pointers
c  for each level are in newstl. also reclaim some space before new
c  allocations.
      call gfixup(lbase, lfnew, nvar, naux)
c
c  merge data structures (newstl and lstart )
c  finish storage allocation, reclaim space, etc. set up boundary
c  flux conservation arrays
c
      do 60 level = lbase, lfine-1
        call prepf(level+1,nvar,naux)
        call prepc(level,nvar)
 60   continue
c
      return
      end
