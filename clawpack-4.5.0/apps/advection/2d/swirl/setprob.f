      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname

      common /compsi/ pi
      common /comvt/ tperiod,pi2

c
c     # compute pi, used in psi.f
      pi = 4.d0 * datan(1.d0)
c
c     # save 2*pi and tperiod in common block for use in b4step2:
c
      pi2 = 2.d0*pi
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      read(7,*) tperiod

      return
      end
