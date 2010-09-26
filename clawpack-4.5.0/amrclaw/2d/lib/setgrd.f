c
c -----------------------------------------------------------------
c
      subroutine setgrd (nvar,cut,naux,dtinit,t0)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical  vtime
      data     vtime/.false./
c     # may as well not bother to calculate time step for error est.
c
c :::::::::::::::::::::::::::: SETGRD :::::::::::::::::::::::::::::::;
c  set up the entire tree/grid structure.  only at this time t = 0
c  can we take advantage of initialization routines.
c  remember that regridding/error estimation needs to have two
c  time steps of soln. values.
c  6/21/05: added dtinit arg. to allow for better choice of initial timestep
c   as discovered by advance/setgrd in first step.
c ::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::::;
c
      dtinit = possk(1)
      if (mxnest .eq. 1) go to 99
c
      levnew =  2
      time   = t0
c
 10   if (levnew .gt. mxnest) go to 30
          levold = levnew - 1
          if (lstart(levold) .eq. 0) go to 30
          lbase  = levold
          lfnew  = lbase
c
c  set up level to be flagged. need a solution t=0,and t=dt.
c  error estimation makes next one at t=2*dt for Richardson.
c
         call advanc(levold,nvar,dtlev,vtime,naux)
         evol = evol + rvol
         rvol = 0.d0
         kfac = 1
         do i = 1, levold-1
           kfac = kfac * kratio(i)
         end do
         dtinit = min(dtinit, dtlev*kfac)
 
c        don't count it in real integration stats
         do 20 level=1,mxnest
 20         rvoll(level) = 0.d0
c
c  flag, cluster, and make new grids
c
         call grdfit(lbase,levold,nvar,naux,cut,time,t0)
         if (newstl(levnew) .ne. 0) lfnew = levnew
c
c  init new level. after each iteration. fix the data structure
c  also reinitalize coarser grids so fine grids can be advanced
c  and interpolate correctly for their bndry vals from coarser grids.
c
         call ginit(newstl(levnew),.true., nvar, naux, t0)
         lstart(levnew) = newstl(levnew)
         lfine = lfnew
         call ginit(lstart(levold),.false., nvar, naux, t0)
c
         levnew = levnew + 1
      go to 10
 30   continue
c
c  switch location of old and new storage for soln. vals, 
c  and reset time to 0.0 (or initial time t0)
c
      if (mxnest .eq. 1) go to 99
c
      lev = 1
 40   if ((lev .eq. mxnest) .or. (lev .gt. lfine))  go to 60
        mptr = lstart(lev)
 50        itemp                = node(store1,mptr)
           node(store1,mptr)    = node(store2,mptr)
           node(store2,mptr)    = itemp
           rnode(timemult,mptr) = t0
           mptr                 = node(levelptr,mptr)
           if (mptr .ne. 0) go to 50
       lev = lev + 1
       go to 40
 60   continue
c
c initial updating so can do conservation check. can do before
c bndry flux arrays set, since don't need them for this
c
      do 65 level = 1, lfine-1
         call update(lfine-level,nvar)
 65   continue
c
c set up boundary flux conservation arrays
c
      do 70 level = 1, lfine-1
         call prepf(level+1,nvar,naux)
         call prepc(level,nvar)
 70   continue
c
 99   continue
c
      return
      end
