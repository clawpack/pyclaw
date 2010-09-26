c
c ------------------------------------------------------------------
c
      subroutine filval(val,mitot,mjtot,hx,hy,lev,time,
     1                  valc,auxc,mic,mjc,
     2                  xleft,xright,ybot,ytop,nvar,
     3                  mptr,ilo,ihi,jlo,jhi,aux,naux,locflip)
 
      implicit double precision (a-h,o-z)

      include "call.i"

      dimension   val(mitot,mjtot,nvar), valc(mic,mjc,nvar)
      dimension   aux(mitot,mjtot,naux), auxc(mic,mjc,naux)
      dimension   dudx(max1d), dudy(max1d)
c
c :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
c
c create and fill coarser (lev-1) patch with one extra coarse cell all
c around, plus the ghost cells . will interpolate from this patch to grid mptr 
c without needing special boundary code. 
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      levc    = lev - 1
      lratiox = intratx(levc)
      lratioy = intraty(levc)
      hxcrse  = hx*lratiox
      hycrse  = hy*lratioy
      xl      = xleft  - hxcrse 
      xr      = xright + hxcrse
      yb      = ybot   - hycrse 
      yt      = ytop   + hycrse
c
c     set integer indices for coarser patch enlarged by 1 cell 
c     (can stick out of domain). proper nesting will insure this one
c     call is sufficient.
      iclo   = ilo/lratiox - 1
      jclo   = jlo/lratioy - 1
      ichi   = (ihi+1)/lratiox - 1 + 1
      jchi   = (jhi+1)/lratioy - 1 + 1
      ng     = 0

c    :::  mcapa  is the capacity function index

      if (mcapa .eq. 0) then   !dont need to copy aux stuff along with soln
        if (xperdom .or. yperdom .or. spheredom) then
          call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,
     &                    locflip)
        else
          call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,1,1)
        endif
      else  ! intersect grids and copy all (soln and aux)
        if (xperdom .or. yperdom .or. spheredom) then
          call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,
     &                  jchi,levc,locflip)
        else
          call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,
     &               levc,1,1)
        endif
      endif
      call bc2amr(valc,auxc,mic,mjc,nvar,naux,
     1            hxcrse,hycrse,levc,time,
     2            xl,xr,yb,yt,
     3            xlower,ylower,xupper,yupper,
     4            xperdom,yperdom,spheredom)
c
c     prepare slopes - use min-mod limiters
c
      do 30 j    = 2, mjc-1
      do 30 ivar = 1, nvar
      do 10 i    = 2, mic-1

         slp = valc(i+1,j,ivar) - valc(i,j,ivar)
         slm = valc(i,j,ivar)   - valc(i-1,j,ivar)
         slopex = dmin1(dabs(slp),dabs(slm))*
     .            dsign(1.0d0,valc(i+1,j,ivar) - valc(i-1,j,ivar))
c        # if there's a sign change, set slope to 0.
         if ( slm*slp .gt. 0.d0) then
           dudx(i) = slopex
         else
           dudx(i) = 0.d0
         endif

         slp = valc(i,j+1,ivar) - valc(i,j,ivar)
         slm = valc(i,j,ivar)   - valc(i,j-1,ivar)
         slopey = dmin1(dabs(slp),dabs(slm))*
     .            dsign(1.0d0,valc(i,j+1,ivar) - valc(i,j-1,ivar))
         if ( slm*slp .gt. 0.d0) then
           dudy(i) = slopey
         else
           dudy(i) = 0.d0
         endif

 10   continue
c
c     interp. from coarse cells to fine grid
c
      do 20 ico = 1,lratiox
      xoff = (float(ico) - .5)/lratiox - .5
         do 20 jco = 1,lratioy
         jfine = (j-2)*lratioy + nghost + jco
         yoff  = (float(jco) - .5)/lratioy - .5
            do 20 i = 2, mic-1
            ifine   = (i-2)*lratiox + nghost + ico
            val(ifine,jfine,ivar) = valc(i,j,ivar)
     1                                + xoff*dudx(i) + yoff*dudy(i)
 20   continue
c
 30   continue

      if (mcapa .ne. 0) then
        call fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,
     &                nvar,naux,levc)
      endif
c
c  overwrite interpolated values with fine grid values, if available.
c
      call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,
     &             jlo-nghost,jhi+nghost,lev,1,1)

      return
      end
