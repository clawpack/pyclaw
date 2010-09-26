      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /param/ gamma, gamma1
c
c     # Set gamma and gamma1 = gamma-1 for Euler equations
c     # Passed to the Riemann solver rp1.f in a common block
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                


      read(7,*) gamma
      gamma1 = gamma - 1.d0

      return
      end
