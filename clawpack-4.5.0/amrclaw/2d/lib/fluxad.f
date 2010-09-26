c
c -------------------------------------------------------
c
      subroutine fluxad(xfluxm,xfluxp,yfluxm,yfluxp,
     1                  svdflx,mptr,mitot,mjtot,
     2                   nvar,lenbc,lratiox,lratioy,ng,dtf,dx,dy)
c

      implicit double precision (a-h,o-z)

      include  "call.i"

c :::::::::::::::::::: FLUXAD ::::::::::::::::::::::::::::::::::
c  save fine grid fluxes  at the border of the grid, for fixing
c  up the adjacent coarse cells. at each edge of the grid, only
c  save the plus or minus fluxes, as necessary. For ex., on
c  left edge of fine grid, it is the minus xfluxes that modify the
c  coarse cell.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      dimension xfluxm(mitot,mjtot,nvar), yfluxm(mitot,mjtot,nvar)
      dimension xfluxp(mitot,mjtot,nvar), yfluxp(mitot,mjtot,nvar)
      dimension svdflx(nvar,lenbc)
 
      nx  = mitot-2*ng
      ny  = mjtot-2*ng
      nyc = ny/lratioy
      nxc = nx/lratiox
 
c ::::: left side saved first
      lind = 0

      do 100 j=1,nyc
         lind = lind + 1
         jfine = (j-1)*lratioy + ng
         do 110 ivar = 1, nvar
            do 120 l=1,lratioy
               svdflx(ivar,lind) = svdflx(ivar,lind) +
     1                             xfluxm(ng+1,jfine+l,ivar)*dtf*dy
c              write(dbugunit,900)lind,xfluxm(1,jfine+l,ivar),
c     .                           xfluxp(1,jfine+l,ivar)
 900           format(' lind ', i4,' m & p ',2e15.7,' svd ',e15.7)
120         continue
110      continue
100   continue
 
c ::::: top side
c      write(dbugunit,*)" saving top side "
      do 200 i=1,nxc
         lind = lind + 1
         ifine = (i-1)*lratiox + ng
         do 210 ivar = 1, nvar
            do 220 l=1,lratiox
               svdflx(ivar,lind) = svdflx(ivar,lind) + 
     1                     yfluxp(ifine+l,mjtot-ng+1,ivar)*dtf*dx
c              write(dbugunit,900)lind,yfluxm(ifine+l,mjtot-ng+1,
c     .                           ivar),yfluxp(ifine+l,mjtot-ng+1,ivar),
c     .                           svdflx(ivar,lind)
220         continue
210      continue
200   continue
 
c ::::: right side
      do 300 j=1,nyc
         lind = lind + 1
         jfine = (j-1)*lratioy + ng
         do 310 ivar = 1, nvar
            do 320 l=1,lratioy
               svdflx(ivar,lind) = svdflx(ivar,lind) + 
     1                     xfluxp(mitot-ng+1,jfine+l,ivar)*dtf*dy
c              write(dbugunit,900)lind,xfluxm(mitot-ng+1,jfine+l,
c                                 ivar),xfluxp(mitot-ng+1,jfine+l,ivar)
320         continue
310      continue
300   continue
 
c ::::: bottom side
c      write(dbugunit,*)" saving bottom side "
      do 400 i=1,nxc
         lind = lind + 1
         ifine = (i-1)*lratiox + ng
         do 410 ivar = 1, nvar
            do 420 l=1,lratiox
               svdflx(ivar,lind) = svdflx(ivar,lind) +
     1                             yfluxm(ifine+l,ng+1,ivar)*dtf*dx
c              write(dbugunit,900)lind,yfluxm(ifine+l,ng+1,ivar),
c     .                      yfluxp(ifine+l,ng+1,ivar),svdflx(ivar,lind)
420         continue
410      continue
400   continue
 
      return
      end
