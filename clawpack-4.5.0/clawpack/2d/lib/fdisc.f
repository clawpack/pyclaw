c
c
c
c     =================================================
      function fdisc(x,y)
c     =================================================

      implicit double precision (a-h,o-z)

c
c     # For computing cell averages for initial data or coefficients that
c     # have a discontinuity along some curve.  

c     # fdisc should be negative to the "left" of the curve and 
c     # positive to the "right".

c     # The cellave routine can then be used to compute the fraction wl of
c     # a grid cell that lies to the "left" of the curve.

c     # Sample code for the case where the curve is the unit circle:

      fdisc = x**2 + y**2 - 1.d0
c
      return
      end
