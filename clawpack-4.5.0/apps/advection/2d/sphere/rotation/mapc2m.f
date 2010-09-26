c
c     ==========================================
      subroutine mapc2m(x1,y1,xp,yp,zp)
c     ==========================================
c
c     Map computational point (x1,y1) to the point (xp,yp,zp)
c     on the sphere
c
c     Points in [-1,1]x[-1,1]  are mapped to the upper hemisphere
c     Points in [-3,-1]x[-1,1] are mapped to the lower hemisphere
c
      implicit real*8(a-h,o-z)

      r1 = 1.d0
      xc = x1
      yc = y1

c     # ghost cell values outside of [-3,1]x[-1,1] get mapped to other
c     # hemisphere:
      if (xc.ge.1.d0) xc = xc-4.d0
      if (xc.lt.-3.d0) xc = xc+4.d0

      if (yc.ge.1.d0) then
         yc = 2.d0 - yc
         xc = -2.d0 - xc
         endif
      if (yc.lt.-1.d0) then
         yc = -2.d0 - yc
         xc = -2.d0 - xc
         endif

      if (xc.lt.-1.d0) then
c         # points in [-3,-1] map to lower hemisphere - reflect about x=-1
c         # to compute x,y mapping and set sgnz appropriately:
          xc = -2.d0 - xc
          sgnz = -1.d0
         else
          sgnz = 1.d0
         endif
      sgnxc = dsign(1.d0,xc)
      sgnyc = dsign(1.d0,yc)

      xc1 = dabs(xc)
      yc1 = dabs(yc)
      d = dmax1(dmax1(xc1,yc1), 1.d-10)

      DD = r1*d*(2.d0 - d) / dsqrt(2.d0)
      R = r1

      center = DD - dsqrt(dmax1(R**2 - DD**2, 0.d0))
      xp = DD/d * xc1
      yp = DD/d * yc1

      if (yc1 .gt. xc1) then
          yp = center + dsqrt(dmax1(R**2 - xp**2, 0.d0))
        else
          xp = center + dsqrt(dmax1(R**2 - yp**2, 0.d0))
        endif

      zp = dsqrt(dmax1(r1**2 - (xp**2 + yp**2), 0.d0))
      xp = xp*sgnxc
      yp = yp*sgnyc
      zp = zp*sgnz
c
      return
      end

