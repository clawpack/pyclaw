      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comrp/ theta
c
c     # Burgers' equation at an angle theta to the grid,
c     #  u_t + cos(theta)*(0.5*u^2)_x + sin(theta)*(0.5*u^2)_y = 0
c
c     # Set theta for angle 
c     # Passed to the Riemann solver rp1.f in a common block
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
      read(7,*) theta

      return
      end

