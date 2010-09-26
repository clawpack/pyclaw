
c
c     =================================================
      double precision function psi(x,y)
c     =================================================
c
c     # stream function 

      implicit double precision (a-h,o-z)
      common /compsi/ pi

      psi = ((dsin(pi*x))**2 * (dsin(pi*y))**2) / pi
c
      return
      end

