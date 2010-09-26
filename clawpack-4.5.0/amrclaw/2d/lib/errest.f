c
c -------------------------------------------------------------
c
      subroutine errest (nvar,naux,lcheck)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

c
c     #   no sense computing new time step if just for error estimation,
c     #   vtime removed from arg list to stepgrid, but new dt not used
 
c :::::::::::::::::::::::::: ERREST :::::::::::::::::::::::::::::::::::
c for all grids at level lcheck:
c  estimate the error by taking a large (2h,2k) step based on the
c  values in the old storage loc., then take one regular (and for
c  now wasted) step based on the new info.   compare using an
c  error relation for a pth order  accurate integration formula.
c  flag error plane as either bad (needs refinement), or good.
c
c  call the regular integrator on a grid twice as coarse.
c  initialize such a grid directly, instead of trick dimensioning.
c  this is to make other l1 type estimates easier to program.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
       hx   = hxposs(lcheck)
       hy   = hyposs(lcheck)
       hx2  = 2.d0*hx
       hy2  = 2.d0*hy
       mptr = lstart(lcheck)
 5     continue
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          locnew = node(store1,mptr)
          locold = node(store2,mptr)
          locaux = node(storeaux,mptr)
          mi2tot = nx/2  + 2*nghost
          mj2tot = ny/2  + 2*nghost
          time   = rnode(timemult,mptr)
          dt     = possk(lcheck)
          tpre   = time - dt
c
c     prepare double the stencil size worth of boundary values,
c            then coarsen them for the giant step integration.
c 
          midub = nx+4*nghost
          mjdub = ny+4*nghost
          locdub = igetsp(midub*mjdub*(nvar+naux))
          locbgc = igetsp(mi2tot*mj2tot*(nvar+naux))
          node(errptr,mptr) = locbgc
          ng2 = 2*nghost

c         # transfer soln. into grid with twice the ghost cells
          call copysol(alloc(locdub),alloc(locold),nvar,mitot,mjtot,
     1              nghost,midub,mjdub,ng2)

c
          if (naux .gt. 0) then
              locaxb = locdub + midub*mjdub*nvar
              locaxc = locbgc + mi2tot*mj2tot*nvar
              xl     = rnode(cornxlo, mptr)
              yb     = rnode(cornylo, mptr)
              maxmx = midub - 4*nghost
              mx = maxmx
              maxmy = mjdub - 4*nghost
              my = maxmy
              call setaux(maxmx,maxmy,2*nghost,mx,my,xl,yb,hx,hy,
     &                    naux,alloc(locaxb))
              call auxcoarsen(alloc(locaxb),midub,mjdub,
     1                     alloc(locaxc),mi2tot,mj2tot,naux,auxtype)
          else
              locaxb = 1
          endif

c         # fill it - use enlarged (before coarsening) aux arrays
          call bound(tpre,nvar,ng2,alloc(locdub),midub,mjdub,mptr,
     1               alloc(locaxb),naux)

c         coarsen by 2 in every direction
          call coarsen(alloc(locdub),midub,mjdub,
     1                 alloc(locbgc),mi2tot,mj2tot,nvar)
          call reclam(locdub,midub*mjdub*(nvar+naux))
c
c
c************************************************************************
c        changed. now locbig filled from spest, always, even if no
c        user estimation routine. reclaim space later tho.   

c We now fill bndry values at time t = time, in preparation
c for calculating the solution on the grid mptr for error estimation.
c 
c          locbig = igetsp(mitot*mjtot*nvar)
c          node(tempptr,mptr) = locbig
c         # straight copy into scratch array so don't mess up latest soln.
c          do 10 i = 1, mitot*mjtot*nvar
c 10          alloc(locbig+i-1) = alloc(locnew+i-1)
c
c          call bound(time,nvar,nghost,alloc(locbig),mitot,mjtot,mptr,
c     1               alloc(locaux),naux)

c************************************************************************
c
       mptr = node(levelptr,mptr)
       if (mptr .ne. 0) go to 5
c
       hx   = hxposs(lcheck)
       hy   = hyposs(lcheck)
       hx2  = 2.d0*hx
       hy2  = 2.d0*hy
       dt   = possk(lcheck)
       dt2    = 2. * dt

       mptr = lstart(lcheck)
 25    continue
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx+ 2*nghost
          mjtot  = ny+ 2*nghost
          mi2tot = nx/2 + 2*nghost
          mj2tot = ny/2 + 2*nghost
c
c         # this scratch storage will be used both for regular and half
c         # sized grid. different dimensions in stepgrid - do not reuse.
          locfp = igetsp(mitot*mjtot*nvar)
          locfm = igetsp(mitot*mjtot*nvar)
          locgp = igetsp(mitot*mjtot*nvar)
          locgm = igetsp(mitot*mjtot*nvar)
          locaux = node(storeaux,mptr)
c
          locbgc = node(errptr,mptr)
c
          locaxc = locbgc+nvar*mi2tot*mj2tot
c should we set to 1 if naux=0?
c
          evol = evol + (nx/2)*(ny/2)
          xlow = rnode(cornxlo,mptr) - nghost*hx2
          ylow = rnode(cornylo,mptr) - nghost*hy2
          call stepgrid(alloc(locbgc),alloc(locfm),alloc(locfp),
     1                alloc(locgm),alloc(locgp),
     2                mi2tot,mj2tot,nghost,
     3                dt2,dtnew2,hx2,hy2,nvar,
     4                xlow,ylow,tpre,mptr,naux,alloc(locaxc))
c
c  the one giant step based on old values is done. now take
c  one regular step based on new values. 
c
      evol   = evol + nx * ny
      xlow   = rnode(cornxlo,mptr) - nghost*hx
      ylow   = rnode(cornylo,mptr) - nghost*hy
      locbig = node(tempptr,mptr)
      loctmp = node(store2,mptr)
c
c **********************************************************************
c ****** changed. spatial error now done from separate routine spest ******
c
c estimate spatial component of error - use old values before
c integration to get accurate boundary gradients
c
c     locerrsp = igetsp(mitot*mjtot)
c     call errsp(alloc(locbig), alloc(locerrsp), mitot,mjtot,
c    &           nvar, mptr,nghost, eprint, outunit)

c **********************************************************************

      call stepgrid(alloc(locbig),alloc(locfm),alloc(locfp),
     1            alloc(locgm),alloc(locgp),
     2            mitot,mjtot,nghost,
     3            dt,dtnew,hx,hy,nvar,
     4            xlow,ylow,time,mptr,naux,alloc(locaux))
c
      call reclam(locfp, mitot*mjtot*nvar)
      call reclam(locfm, mitot*mjtot*nvar)
      call reclam(locgp, mitot*mjtot*nvar)
      call reclam(locgm, mitot*mjtot*nvar)
c
      
      call errf1(alloc(locbig),nvar,
     1          alloc(locbgc),mptr,mi2tot,mj2tot,
     2          mitot,mjtot,alloc(loctmp))

      call reclam(locbgc,mi2tot*mj2tot*(nvar+naux))

c      reclaiming locbig is now done from bufnst, since this routine 
c      will not be called if tol<0
c      call reclam(locbig,mitot*mjtot*nvar)
c
      mptr = node(levelptr, mptr)
      if (mptr .ne. 0) go to 25
c
      return
      end
