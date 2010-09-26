      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comrp/ u
      common /comic/ beta
c
c     # Set the velocity for scalar advection
c     # This value is passed to the Riemann solver rp1.f in a common block
c
c     # Set the width of the initial Gaussian pulse
c     # beta is passed to qinit.f in comic
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
      read(7,*) u
      read(7,*) beta

      return
      end
