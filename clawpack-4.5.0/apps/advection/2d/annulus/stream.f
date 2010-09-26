c
c     =================================================
      function stream(x,y)
c     =================================================
c
c     # Stream function in physical space (x,y).
c     # Clockwise rotation, rotates fully in time 1.

      implicit double precision (a-h,o-z)

      stream = 3.1415926535897931d0 *(x**2 + y**2)
c
      return
      end
