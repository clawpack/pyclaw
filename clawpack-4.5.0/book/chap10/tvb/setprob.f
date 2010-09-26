      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comrp/ u
      common /comic/ beta, freq
      common /comlim/ phiM,phiMdx2
c
c     # Set the velocity for scalar advection
c     # This value is passed to the Riemann solver rp1.f in a common block
c
c     # Set the width of the Gaussian for the wave packet
c     # beta is passed to qinit.f in comic
c
c     # Set the frequency of the wave packet
c     # freq is passed to qinit.f in comic
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
      read(7,*) u
      read(7,*) beta
      read(7,*) freq
      read(7,*) phiM

      return
      end
