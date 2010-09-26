c
c     =====================================================
      subroutine mapc2p(xc,yc,xp,yp)
c     =====================================================
c
c     # on input,  (xc,yc) is a computational grid point
c     # on output, (xp,yp) is corresponding point in physical space
c
      implicit double precision (a-h,o-z)
c
c     # Polar coordinates, xc = r,  yc = theta
      xp = xc * dcos(yc) 
      yp = xc * dsin(yc)
      return
      end
