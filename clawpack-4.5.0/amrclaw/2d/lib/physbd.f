c
c ------------------------------------------------------------------
c
      subroutine physbd(val,aux,nrow,ncol,nvar,naux,
     1                  hx, hy, level, time,
     2                  xleft, xright, ybot, ytop,
     3                  xlower,ylower,xupper,yupper,
     4                  xperiodic,    yperiodic)

 
c
c
c :::::::::: PHYSBD ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     Take a grid patch with mesh widths hx,hy, of dimensions nrow by
c     ncol,  and set the values of any piece of
c     of the patch which extends outside the physical domain 
c     using the boundary conditions. 
c
c     The corners of the grid patch are at 
c        (xleft,ybot)  --  lower left corner
c        (xright,ytop) --  upper right corner
c
c     The physical domain itself is a rectangle bounded by
c        (xlower,ylower)  -- lower left corner
c        (xupper,yupper)  -- upper right corner
c     
c     the picture is the following: 
c
c               _____________________ (xupper,yupper)
c              |                     |  
c          _________ (xright,ytop)   |
c          |   |    |                |
c          |   |    |                |
c          |   |    |                |
c          |___|____|                |
c (xleft,ybot) |                     |
c              |                     |
c              |_____________________|
c   (xlower,ylower)
c        
c
c     Any cells that lie outside the physical domain are ghost cells whose
c     values should be set in this routine.  This is tested for by comparing
c     xleft with xlower to see if values need to be set at the left, as in
c     the figure above, and similarly at the other boundaries.
c
c     Patches are guaranteed to have at least 1 row of cells filled
c     with interior values so it is possible to  extrapolate. 
c     Fix trimbd if you want more than 1 row pre-set.
c
c     Make sure the order the boundaries are specified is correct
c     so that diagonal corner cells are also properly taken care of.
c
c     Periodic boundaries are set before calling this routine, so you
c     can safely extrapolate there.  Don't overwrite them!
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      implicit double precision (a-h,o-z)
      logical xperiodic, yperiodic

      include  "cuser.i"

      dimension val(nrow,ncol,nvar), aux(nrow,ncol,naux)


c     right boundary 
      if (xright .gt. xupper+hxmarg) then
        nxr  = (xright - xupper + hxmarg)/hx
        nxr  = nxr - 1
        ibeg = max0(nrow-nxr, 1)
             do 11 i    = ibeg, nrow
             do 11 j    = 1, ncol
             do 11 ivar = 1, nvar
               val(i,j,ivar) = val(2*ibeg-1-i,j,ivar) 
11           continue
      endif

c     bottom boundary - reflecting (for ramp part, else inflow
      if (ybot .lt. ylower-hymarg) then
        nyb = (ylower+hymarg-ybot)/hy
         do 12 i=1,nrow
         do 12 j=1,nyb     
            val(i,j,ivar)
  12     continue
c
      endif
 

c     top boundary - reflecting
      if (ytop .gt. yupper+hymarg) then
          nyt = (ytop - yupper + hymarg)/hy
          jbeg = max0(ncol-nyt+1, 1)
          x0      =  disp + 10.d0*time/cos(pi/6.d0)
          y0      =  ylower

           do 13 j=jbeg,ncol
             ycen   =  ybot + float(j-.5d0)*hy

              do 13 i=1,nrow
                xcen   = xleft + float(i-.5d0)*hx

                call cellave(xcen-hx/2.d0,ycen-hy/2.d0,hx,hy,wl)

                rho = (1.d0-wl)*rhoamb + wl*rhoshk
                u   = (1.d0-wl)*uamb + wl*ushk
                v   = (1.d0-wl)*vamb + wl*vshk
                p   = (1.d0-wl)*pamb + wl*pshk

                val(i,j,1) =  rho
                val(i,j,2) =  rho * u
                val(i,j,3) =  rho * v
                val(i,j,4) =  p/gamma1 + .5d0*rho*(u*u+v*v)

  13       continue

      endif

c     left boundary - inflow
      if (xleft .lt. xlower-hxmarg) then
        nxl = (xlower+hxmarg-xleft)/hx
            do 10 i = 1, nxl
            do 10 j = 1, ncol
              val(i,j,1) = rhoshk 
              val(i,j,2) = rhoshk * ushk
              val(i,j,3) = rhoshk * vshk
              val(i,j,4) = rhoshk*(eshk+.5d0*(ushk**2+vshk**2))
 10         continue
      endif

      return
      end
      return
      end
