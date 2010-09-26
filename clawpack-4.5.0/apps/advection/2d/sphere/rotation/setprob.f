      subroutine setprob
      implicit double precision (a-h,o-z)

      common /comaxis/ xaxis,yaxis,zaxis

c     # axis of rotation and speed:

      xaxis = 0.d0
      yaxis = 0.d0
      zaxis = 6.2832d0

      close(unit=7)

      return
      end
