      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comrp/ grav
      common /comsrc/ ndim
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /comic/ hin,hout
c
c     # Set the material parameters for the acoustic equations
c     # Passed to the Riemann solver rp1.f in a common block
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
c     # ndim = space dimensions (2 = cylindrical symmetry, 3 = spherical)
      read(7,*) ndim

c     # gravitational constant:
      read(7,*) grav

c     # data for radial dam-break problem:
      read(7,*) x0,y0,r0
      read(7,*) hin,hout
c
      return
      end
