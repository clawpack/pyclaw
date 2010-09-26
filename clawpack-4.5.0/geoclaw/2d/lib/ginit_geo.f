c
c -------------------------------------------------------------
c
      subroutine ginit(msave, first, nvar, naux, t0)
c
      implicit double precision (a-h,o-z)
      logical first

      include  "call.i"

c ::::::::::::::::::::::::::::: GINIT ::::::::::::::::::::::::
c
c  initializes soln on all grids at 'level'  by calling qinit
c  if first = true, (first call to init), then allocate the
c  soln storage area too, else was already allocated.
c
c :::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::

      if (msave .eq. 0) go to 99
      level = node(nestlevel,msave)
      hx    = hxposs(level)
      hy    = hyposs(level)
      mptr  = msave

 10       nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot   = nx + 2*nghost
          mjtot   = ny + 2*nghost
          corn1   = rnode(cornxlo,mptr)
          corn2   = rnode(cornylo,mptr)
          if(.not. (first)) go to 20
              loc                 = igetsp(mitot*mjtot*nvar)
              node(store1,mptr)   = loc
              if (naux .gt. 0) then
                locaux              = igetsp(mitot*mjtot*naux)
                maxmx = mitot - 2*nghost
                mx = maxmx
                maxmy = mjtot - 2*nghost
                my = maxmy
                call setaux(maxmx,maxmy,nghost,mx,my,corn1,corn2,hx,hy,
     &                    naux,alloc(locaux))
              else 
                locaux = 1
              endif
              node(storeaux,mptr) = locaux
              if (level .lt. mxnest) then
                loc2              = igetsp(mitot*mjtot*nvar)
                node(store2,mptr) = loc2
              endif
              rnode(timemult, mptr) = t0
              go to 30
 20       continue
c
c  if 2nd time through, put initial values in store2 so finer grids
c  can be advanced with interpolation of their boundary values.
c  new time soln should still be in location store1.
c
          loc     = node(store2,mptr)
          locaux  = node(storeaux,mptr)
c
   30     continue
          maxmx = mitot - 2*nghost
          mx = maxmx
          maxmy = mjtot - 2*nghost
          my = maxmy
          call qinit(maxmx,maxmy,nvar,nghost,mx,my,corn1,corn2,hx,hy,
     &               alloc(loc),naux,alloc(locaux))


c         call qinit(alloc(loc),nvar,mitot,mjtot,
c    1               corn1-nghost*hx,corn2-nghost*hy,hx,hy,
c    2               alloc(locaux),naux)
c
          mptr  = node(levelptr, mptr)
      if (mptr .ne. 0) go to 10
c
c
 99   return
      end
c
c ----------------------------------------------------------------------------------------------------------
c
      double precision function getMaxSpeed(val,mitot,mjtot,nvar,
     &                          aux,naux,nghost,hx,hy)

      use geoclaw_module
      implicit double precision (a-h,o-z)
      dimension   val(mitot,mjtot,nvar), aux(mitot,mjtot,naux)


      sp_over_h = 0.d0   ! compute max speed over h, since dx may not equal dy
      if (icoordsys .eq. 2) then
        do j = nghost+1, mjtot-nghost
         ymetric = Rearth*pi/180.d0
         hyphys = ymetric*hy

         do i = nghost+1, mitot-nghost
            xmetric = cos(aux(i,j,3))*Rearth*pi/180.d0
            hxphys = xmetric * hx
            h  = val(i,j,1)
            if (h .gt. drytolerance) then
              u  = val(i,j,2)/h
              v  = val(i,j,3)/h
            else
              u = 0.d0
              v = 0.d0
            endif
            sig = sqrt(grav*h)
            sp_over_h = max((abs(u)+sig)/hxphys,(abs(v)+sig)/hyphys,
     &                      sp_over_h)
         end do
         end do
      else  ! speeds in cartesian coords, no metrics needed
         do j = nghost+1, mjtot-nghost
         do i = nghost+1, mitot-nghost
            h  = val(i,j,1)
            if (h .gt. drytolerance) then
              u  = val(i,j,2)/h
              v  = val(i,j,3)/h
            else
              u = 0.d0
              v = 0.d0
            endif
            sig = sqrt(grav*h)
            sp_over_h = max((abs(u)+sig)/hx,(abs(v)+sig)/hy,sp_over_h)
         end do
         end do
      endif

       getMaxSpeed = sp_over_h

       return
       end
