      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comrp/ u
      common /comsrc/ dcoef
c
c     # Set the velocity for scalar advection
c     # This value is passed to the Riemann solver rp1.f in a common block
c
c     # Set the diffusion coefficient for the source term
c     # dcoef is passed to src1.f in comsrc
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
      read(7,*) u
      read(7,*) dcoef

      return
      end
