c
c --------------------------------------------------------------
c
      subroutine bound(time,nvar,ng,valbig,mitot,mjtot,mptr,
     1                 aux,naux)

c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension valbig(mitot,mjtot,nvar), aux(mitot,mjtot,naux)

c
c  :::::::::::::: BOUND :::::::::::::::::::::::::::::::::::::::::::
c     This routine sets the boundary values for a given grid 
c     at level level.
c     We are setting the values for a strip ng zones wide all
c     the way around the border, in 4 rectangular strips.
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted
c     directly into the enlarged valbig array.
c
c     This routine calls the routine filpatch
c     which for any block of mesh points on a given level,
c     intersects that block with all grids on that level and with
c     the physical boundaries, copies the values into the
c     appropriate intersecting regions, and interpolates the remaining
c     cells from coarser grids as required.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      xleft  = rnode(cornxlo, mptr)
      xright = rnode(cornxhi, mptr)
      ybot   = rnode(cornylo, mptr)
      ytop   = rnode(cornyhi, mptr)
      ilo    = node(ndilo, mptr)
      ihi    = node(ndihi, mptr)
      jlo    = node(ndjlo, mptr)
      jhi    = node(ndjhi, mptr)
      level  = node(nestlevel, mptr)
      hx     = hxposs(level)
      hy     = hyposs(level)


c     left boundary

      xl = xleft - ng*hx
      xr = xleft
      yb = ybot 
      yt = ytop


      if ((xl .lt. xlower) .and. xperdom) then
        call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                    1,ng+1,
     2                    ilo-ng,ilo-1,jlo,jhi)
      else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                1,ng+1,
     2                ilo-ng,ilo-1,jlo,jhi)
      endif


c     right boundary

      xl = xright
      xr = xright + ng*hx
      yb = ybot
      yt = ytop

      if ((xr .gt. xupper) .and. xperdom) then
        call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                    mitot-ng+1,ng+1,
     2                    ihi+1,ihi+ng,jlo,jhi)
      else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                mitot-ng+1,ng+1,
     2                ihi+1,ihi+ng,jlo,jhi)
      endif


c       bottom boundary
        xl = xleft  - ng*hx
        xr = xright + ng*hx
        yb = ybot - ng*hy
        yt = ybot
        if (((yb .lt. ylower) .and. (yperdom .or. spheredom)) .or. 
     1   (((xl .lt. xlower) .or. (xr .gt. xupper)) .and. xperdom) ) then
           call prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                      1,1,
     2                      ilo-ng,ihi+ng,jlo-ng,jlo-1)
        else
           call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                   1,1,
     2                   ilo-ng,ihi+ng,jlo-ng,jlo-1)
        endif

c       top boundary
        xl = xleft - ng*hx
        xr = xright + ng*hx
        yb = ytop
        yt = ytop + ng*hy
        if (((yt .gt. yupper) .and. (yperdom .or. spheredom)) .or. 
     1   (((xl .lt. xlower) .or. (xr .gt. xupper)) .and. xperdom) ) then
          call prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                     1,mjtot-ng+1,
     2                     ilo-ng,ihi+ng,jhi+1,jhi+ng)
        else
          call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                  1,mjtot-ng+1,
     2                  ilo-ng,ihi+ng,jhi+1,jhi+ng)
        endif

c
c external boundary conditions   THIS IS NOW DONE IN THE FILPATCHES
c (in the recursive filrecur.f, since filpatches had to call bc2amr,
c  have them all do it).
c

!--      if (.not. (xperdom .and. (yperdom .or. spheredom) )) then
!--        xl = xleft  - ng*hx
!--        yb = ybot   - ng*hy
!--        xr = xright + ng*hx
!--        yt = ytop   + ng*hy
!--
!--        call bc2amr( valbig,aux,mitot,mjtot,nvar,naux,
!--     1               hx,hy,level,time,xl,xr,yb,yt,
!--     3               xlower,ylower,xupper,yupper,
!--     4               xperdom,yperdom,spheredom)
!--      endif
c
      return
      end
